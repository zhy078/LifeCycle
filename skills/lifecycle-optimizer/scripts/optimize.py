#!/usr/bin/env python3
import argparse
import itertools
import json
import math
import random
import re
import shutil
import subprocess
import sys
import time
from pathlib import Path


PARAM_LINE_RE = re.compile(r"^(?P<name>\w+)\s*=\s*[^;]+;\s*$")


def resolve_output_dir(raw_path: str) -> Path:
    p = Path(raw_path)
    if p.is_absolute():
        return p
    repo_root = Path(__file__).resolve().parents[3]
    return (repo_root / p).resolve()


def repo_root() -> Path:
    return Path(__file__).resolve().parents[3]


def load_config(path: Path):
    with path.open("r", encoding="utf-8") as f:
        return json.load(f)


def validate(cfg):
    if "constraints" not in cfg or "search_space" not in cfg:
        raise ValueError("config must include constraints and search_space")
    c = cfg["constraints"]
    required = ["rho", "delta", "psi", "mu", "sigr"]
    for k in required:
        if k not in cfg["search_space"] or not cfg["search_space"][k]:
            raise ValueError(f"missing search_space.{k}")

    def in_range(name, v, lo, hi):
        if not (lo <= v <= hi):
            raise ValueError(f"{name} out of range: {v}. expected [{lo}, {hi}]")

    for rho in cfg["search_space"]["rho"]:
        in_range("rho", rho, c.get("rho_min", 1.0), c.get("rho_max", 1000.0))
    for d in cfg["search_space"]["delta"]:
        in_range("delta", d, c.get("delta_min", 0.0), c.get("delta_max", 1.0))
    for psi in cfg["search_space"]["psi"]:
        in_range("psi", psi, c.get("psi_min", 0.0), c.get("psi_max", 1000.0))
    for mu in cfg["search_space"]["mu"]:
        in_range("mu", mu, c.get("mu_min", -1.0), c.get("mu_max", 1.0))
    for sigr in cfg["search_space"]["sigr"]:
        in_range("sigr", sigr, c.get("sigr_min", 0.0), c.get("sigr_max", 10.0))


def score_proxy(params, objective):
    rho = params["rho"]
    delta = params["delta"]
    psi = params["psi"]
    mu = params["mu"]
    sigr = params["sigr"]
    baseline = 10 * math.log1p(max(mu, -0.99) + 1.0) + 50 * delta - 0.7 * rho - 5 * sigr + 2 * psi
    if objective == "maximize_terminal_wealth":
        return baseline
    if objective == "maximize_risk_adjusted":
        return baseline - 3 * sigr
    if objective == "minimize_volatility_penalty":
        return -sigr + 0.2 * mu + 0.1 * psi
    return baseline


def octave_available():
    return shutil.which("octave") is not None


def patch_lifecycle_script(src: Path, dst: Path, params: dict, fixed: dict, fast_mode: bool):
    replacements = {
        "rho": params["rho"],
        "delta": params["delta"],
        "psi": params["psi"],
        "mu": params["mu"],
        "sigr": params["sigr"],
        "r": fixed.get("r", 1.015),
        "tb": fixed.get("tb", 20),
        "tr": fixed.get("tr", 66),
        "td": fixed.get("td", 100),
        "nsim": fixed.get("nsim", 10000),
    }
    if fast_mode:
        replacements.update({
            "na": fixed.get("na", 10),
            "ncash": fixed.get("ncash", 10),
            "nsim": fixed.get("nsim", 1000),
        })

    lines = src.read_text(encoding="utf-8").splitlines()
    out_lines = []
    for ln in lines:
        m = PARAM_LINE_RE.match(ln.strip())
        if m and m.group("name") in replacements:
            name = m.group("name")
            val = replacements[name]
            if isinstance(val, str):
                out_lines.append(f"{name} = {val};")
            else:
                out_lines.append(f"{name} = {val};")
        else:
            out_lines.append(ln)
    dst.write_text("\n".join(out_lines) + "\n", encoding="utf-8")


def parse_terminal_metric(run_dir: Path):
    cwy = run_dir / "CWY.txt"
    if cwy.exists():
        lines = [ln.strip() for ln in cwy.read_text(encoding="utf-8", errors="ignore").splitlines() if ln.strip()]
        if lines:
            parts = lines[-1].replace(",", " ").split()
            if len(parts) >= 2:
                try:
                    return float(parts[1])
                except ValueError:
                    pass

    year80 = run_dir / "year80.txt"
    if not year80.exists():
        return None
    lines = [ln.strip() for ln in year80.read_text(encoding="utf-8", errors="ignore").splitlines() if ln.strip()]
    if not lines:
        return None
    parts = lines[-1].replace(",", " ").split()
    if not parts:
        return None
    try:
        return float(parts[0])
    except ValueError:
        return None


def _read_numeric_pairs(path: Path):
    pairs = []
    for raw in path.read_text(encoding="utf-8", errors="ignore").splitlines():
        line = raw.strip()
        if not line:
            continue
        parts = line.replace(",", " ").split()
        if len(parts) < 2:
            continue
        try:
            pairs.append((float(parts[0]), float(parts[1])))
        except ValueError:
            continue
    return pairs


def parse_year_policy(path: Path):
    rows = _read_numeric_pairs(path)
    if not rows or len(rows) % 3 != 0:
        return None
    block = len(rows) // 3

    def build_segment(segment_rows):
        return [
            {"value": value, "cash": cash}
            for value, cash in segment_rows
        ]

    return {
        "alpha": build_segment(rows[:block]),
        "consumption": build_segment(rows[block:2 * block]),
        "value": build_segment(rows[2 * block:]),
    }


def summarize_policy_curve(entries):
    if not entries:
        return None
    mid = len(entries) // 2
    selected = {
        "low_wealth": entries[0],
        "mid_wealth": entries[mid],
        "high_wealth": entries[-1],
    }
    values = [row["value"] for row in entries]
    return {
        "points": selected,
        "min_value": min(values),
        "max_value": max(values),
    }


def build_policy_summary(run_dir: Path, fixed: dict):
    year_files = sorted(run_dir.glob("year*.txt"))
    if not year_files:
        return None

    tb = int(fixed.get("tb", 20))
    tr = int(fixed.get("tr", 66))
    td = int(fixed.get("td", 100))
    per_age = []
    for idx, year_file in enumerate(year_files):
        parsed = parse_year_policy(year_file)
        if not parsed:
            continue
        age = min(tb + idx, td - 1)
        per_age.append({
            "age": age,
            "phase": "working" if age < tr else "retired",
            "alpha": summarize_policy_curve(parsed["alpha"]),
            "consumption": summarize_policy_curve(parsed["consumption"]),
            "value": summarize_policy_curve(parsed["value"]),
        })

    if not per_age:
        return None

    pivot_ages = []
    for desired in [tb, min(tr - 1, per_age[-1]["age"]), tr, per_age[-1]["age"]]:
        match = next((row for row in per_age if row["age"] == desired), None)
        if match and all(existing["age"] != match["age"] for existing in pivot_ages):
            pivot_ages.append(match)

    avg_alpha_mid = sum(row["alpha"]["points"]["mid_wealth"]["value"] for row in per_age) / len(per_age)
    avg_cons_mid = sum(row["consumption"]["points"]["mid_wealth"]["value"] for row in per_age) / len(per_age)

    return {
        "age_span": {"start": per_age[0]["age"], "end": per_age[-1]["age"]},
        "pivot_ages": pivot_ages,
        "average_mid_wealth_alpha": avg_alpha_mid,
        "average_mid_wealth_consumption": avg_cons_mid,
        "per_age": per_age,
    }


def render_policy_markdown(policy_summary):
    if not policy_summary:
        return ""

    lines = [
        "## Policy Path Summary",
        "",
        f"- age span: {policy_summary['age_span']['start']} to {policy_summary['age_span']['end']}",
        f"- average mid-wealth risky share: {policy_summary['average_mid_wealth_alpha']:.4f}",
        f"- average mid-wealth consumption: {policy_summary['average_mid_wealth_consumption']:.4f}",
        "",
        "### Pivot Ages",
        "",
    ]

    for row in policy_summary["pivot_ages"]:
        alpha = row["alpha"]["points"]
        cons = row["consumption"]["points"]
        lines.append(
            f"- age {row['age']} ({row['phase']}): "
            f"alpha low/mid/high = {alpha['low_wealth']['value']:.4f}/{alpha['mid_wealth']['value']:.4f}/{alpha['high_wealth']['value']:.4f}; "
            f"consumption low/mid/high = {cons['low_wealth']['value']:.4f}/{cons['mid_wealth']['value']:.4f}/{cons['high_wealth']['value']:.4f}"
        )

    return "\n".join(lines) + "\n"


def build_lifecycle_checkpoints(policy_summary):
    if not policy_summary:
        return None

    checkpoints = []
    for row in policy_summary.get("pivot_ages", []):
        wealth_bands = {}
        for band in ["low_wealth", "mid_wealth", "high_wealth"]:
            wealth_bands[band] = {
                "alpha": row["alpha"]["points"][band]["value"],
                "consumption": row["consumption"]["points"][band]["value"],
                "cash": row["alpha"]["points"][band]["cash"],
            }
        checkpoints.append({
            "age": row["age"],
            "phase": row["phase"],
            "mid_wealth_alpha": wealth_bands["mid_wealth"]["alpha"],
            "mid_wealth_consumption": wealth_bands["mid_wealth"]["consumption"],
            "wealth_bands": wealth_bands,
        })

    if not checkpoints:
        return None

    return {
        "age_span": policy_summary["age_span"],
        "checkpoints": checkpoints,
    }


def describe_risk_level(alpha_value):
    if alpha_value >= 0.6:
        return "高"
    if alpha_value >= 0.3:
        return "中"
    return "低"


def phase_label(age, retirement_age):
    if age < 35:
        return "早期积累期"
    if age < retirement_age:
        return "临近退休准备期" if retirement_age - age <= 5 else "中期积累期"
    if age < 80:
        return "退休提款期"
    return "高龄保障期"


def portfolio_mix_from_alpha(alpha_value):
    equity = int(round(alpha_value * 100))
    stabilizer = max(0, 100 - equity)
    if equity >= 60:
        bucket = "成长型组合"
    elif equity >= 30:
        bucket = "平衡型组合"
    else:
        bucket = "稳健型组合"
    return {
        "bucket": bucket,
        "equity_pct": equity,
        "stabilizer_pct": stabilizer,
    }


def portfolio_action(stage, phase):
    if phase == "working" and stage == "早期积累期":
        return "以长期增值为主，权益类资产可以相对更高，持续积累长期资金。"
    if phase == "working":
        return "逐步降低组合波动，为退休前后的提款需求预留更稳定的资产。"
    if stage == "退休提款期":
        return "从追求收益转向现金流稳定，优先保证提款连续性和回撤可控。"
    return "进一步降低风险资产波动暴露，优先保留流动性和保障性资产。"


def wealth_band_label(band):
    return {
        "low_wealth": "低财富",
        "mid_wealth": "中财富",
        "high_wealth": "高财富",
    }.get(band, band)


def choose_wealth_band(wealth_bands, current_wealth):
    if current_wealth is None or not wealth_bands:
        return "mid_wealth"
    return min(
        wealth_bands.keys(),
        key=lambda name: abs(float(wealth_bands[name]["cash"]) - float(current_wealth)),
    )


def build_advisor_insights(best, fixed, client_profile=None):
    checkpoints = (best.get("lifecycle_checkpoints") or {}).get("checkpoints") or []
    if not checkpoints:
        return None

    client_profile = client_profile or {}
    current_wealth = client_profile.get("current_wealth")
    retirement_age = int(fixed.get("tr", 65))
    lines = []
    for row in checkpoints:
        wealth_bands = row.get("wealth_bands") or {}
        selected_band = choose_wealth_band(wealth_bands, current_wealth)
        selected = wealth_bands.get(selected_band, {
            "alpha": row["mid_wealth_alpha"],
            "consumption": row["mid_wealth_consumption"],
            "cash": None,
        })
        risk_level = describe_risk_level(selected["alpha"])
        stage = phase_label(row["age"], retirement_age)
        portfolio = portfolio_mix_from_alpha(selected["alpha"])
        if row["phase"] == "working":
            client_message = (
                f"客户处于{stage}，按照{wealth_band_label(selected_band)}路径，建议以{portfolio['bucket']}为主，当前重点是长期增值与退休前准备。"
            )
        else:
            client_message = (
                f"客户已进入{stage}，按照{wealth_band_label(selected_band)}路径，建议以{portfolio['bucket']}承接退休后的现金流需求，减少大幅回撤风险。"
            )
        lines.append({
            "age": row["age"],
            "phase": row["phase"],
            "stage": stage,
            "risk_level": risk_level,
            "portfolio": portfolio,
            "selected_wealth_band": selected_band,
            "selected_band_label": wealth_band_label(selected_band),
            "selected_cash_reference": selected.get("cash"),
            "selected_alpha": selected["alpha"],
            "selected_consumption": selected["consumption"],
            "wealth_bands": wealth_bands,
            "mid_wealth_alpha": row["mid_wealth_alpha"],
            "mid_wealth_consumption": row["mid_wealth_consumption"],
            "portfolio_action": portfolio_action(stage, row["phase"]),
            "client_message": client_message,
        })
    return lines


def render_customer_manager_report(best, report_meta, fixed, client_profile=None):
    checkpoints = (best.get("lifecycle_checkpoints") or {}).get("checkpoints") or []
    if not checkpoints:
        return "# Customer Manager Report\n\n暂无可用于客户沟通的生命周期检查点。\n"

    client_profile = client_profile or {}
    current_wealth = client_profile.get("current_wealth")
    insights = build_advisor_insights(best, fixed, client_profile) or []
    retirement_age = int(fixed.get("tr", 65))
    death_age = int(fixed.get("td", 100))
    objective = best.get("objective", report_meta.get("objective", "maximize_terminal_wealth"))
    params = best.get("params", {})

    first = insights[0]
    last = insights[-1]
    summary = (
        f"客户当前处于`{first['stage']}`，按照`{first['selected_band_label']}`路径，建议以`{first['portfolio']['bucket']}`作为主配置思路，"
        f"并在`{retirement_age}岁`退休后逐步转向更稳健的提款型组合。"
    )

    lines = [
        "# Customer Manager Report",
        "",
        "## 客户结论",
        "",
        summary,
        "",
        f"- 生命周期区间: {fixed.get('tb', 20)}岁开始工作，{retirement_age}岁退休，{death_age}岁寿命终点",
        f"- 模型目标: {objective}",
        f"- 最优场景分数: {best.get('score', 0.0):.6f}",
        f"- 财富输入: {'未提供，当前报告默认展示所选财富路径并保留其他财富档位供比对' if current_wealth is None else f'current_wealth={current_wealth}，已按最接近的 policy function 节点取值'}",
        f"- 建议起始组合: {first['portfolio']['bucket']}，权益类约{first['portfolio']['equity_pct']}%，稳健类约{first['portfolio']['stabilizer_pct']}%",
        f"- 模型参数附录: rho={params.get('rho')}, delta={params.get('delta')}, psi={params.get('psi')}, mu={params.get('mu')}, sigr={params.get('sigr')}",
        "",
        "## 推荐组合",
        "",
    ]

    for row in insights:
        lines.append(
            f"- {row['age']}岁 | {row['stage']} | 推荐 `{row['portfolio']['bucket']}` | "
            f"权益类约{row['portfolio']['equity_pct']}%，稳健类约{row['portfolio']['stabilizer_pct']}% | "
            f"财富路径: {row['selected_band_label']} (参考 cash={row['selected_cash_reference']}) | "
            f"原因: {row['portfolio_action']}"
        )

    lines.extend([
        "",
        "## 财富敏感度",
        "",
    ])

    for row in insights:
        wealth_bands = row["wealth_bands"]
        lines.append(
            f"- {row['age']}岁: "
            f"低财富 alpha={wealth_bands['low_wealth']['alpha']:.4f}, "
            f"中财富 alpha={wealth_bands['mid_wealth']['alpha']:.4f}, "
            f"高财富 alpha={wealth_bands['high_wealth']['alpha']:.4f}"
        )

    lines.extend([
        "",
        "## 调仓方向",
        "",
    ])

    for current, nxt in zip(insights, insights[1:]):
        lines.append(
            f"- 从{current['age']}岁到{nxt['age']}岁: 组合从 `{current['portfolio']['bucket']}` 过渡到 "
            f"`{nxt['portfolio']['bucket']}`，权益类建议由约{current['portfolio']['equity_pct']}% 调整为约{nxt['portfolio']['equity_pct']}%。"
        )

    lines.extend([
        "",
        "## 客户经理沟通话术",
        "",
    ])

    for row in insights:
        lines.append(f"- {row['age']}岁节点: {row['client_message']}")

    lines.extend([
        "",
        "## 风险提示",
        "",
        f"- 退休年龄是组合切换的关键拐点。如果客户真实退休计划不是 {retirement_age} 岁，建议重新测算。",
        "- 本报告给的是配置方向，不是具体基金或产品名单；落地时还要结合客户风险评级和产品准入。",
        "- 如果客户近期有大额支出、医疗安排或提前退休计划，客户经理应人工覆盖模型建议。",
        "",
        "## 下一步动作",
        "",
        "- 先确认客户真实退休年龄、退休后年支出和可接受回撤区间。",
        "- 再把客户当前持仓映射成成长型、平衡型、稳健型三类资产，看和建议组合差多少。",
        "- 最后形成可执行调仓单，明确哪些资产保留、哪些逐步降低、哪些用于退休现金流准备。",
        "",
        "## 关键观察",
        "",
        f"- 退休前主线: {insights[0]['age']}岁到{insights[1]['age']}岁，组合核心仍是 `{insights[0]['portfolio']['bucket']}`。",
        f"- 退休切换: {retirement_age - 1}岁到{retirement_age}岁，建议开始把沟通重点从增值转到提款稳定。",
        f"- 长寿阶段: {last['age']}岁时，权益类建议降到约{last['portfolio']['equity_pct']}%，更强调保障属性。",
        "",
        f"_生成说明: 共评估 {report_meta.get('total_scenarios')} 个场景，模型运行 {report_meta.get('elapsed_seconds')} 秒。_",
        "",
    ])
    return "\n".join(lines)


def describe_risk_level(alpha_value):
    if alpha_value >= 0.6:
        return "high"
    if alpha_value >= 0.3:
        return "medium"
    return "low"


def phase_label(age, retirement_age):
    if age < 35:
        return "early_accumulation"
    if age < retirement_age:
        return "pre_retirement" if retirement_age - age <= 5 else "mid_accumulation"
    if age < 80:
        return "retirement_income"
    return "longevity_protection"


def phase_label_cn(phase_key):
    return {
        "early_accumulation": "早期积累期",
        "mid_accumulation": "中期积累期",
        "pre_retirement": "临近退休准备期",
        "retirement_income": "退休提款期",
        "longevity_protection": "高龄保障期",
    }[phase_key]


def risk_label_cn(risk_level):
    return {
        "high": "高",
        "medium": "中",
        "low": "低",
    }[risk_level]


def portfolio_mix_from_alpha(alpha_value):
    equity = int(round(alpha_value * 100))
    stabilizer = max(0, 100 - equity)
    if equity >= 60:
        bucket = "growth"
    elif equity >= 30:
        bucket = "balanced"
    else:
        bucket = "conservative"
    return {
        "bucket": bucket,
        "equity_pct": equity,
        "stabilizer_pct": stabilizer,
    }


def bucket_label_cn(bucket):
    return {
        "growth": "成长型组合",
        "balanced": "平衡型组合",
        "conservative": "稳健型组合",
    }[bucket]


def portfolio_action(stage, phase):
    if phase == "working" and stage == "early_accumulation":
        return "以长期增值为主，权益类资产可以相对更高，持续积累长期资金。"
    if phase == "working":
        return "逐步降低组合波动，为未来阶段的资金安排预留更稳定的资产。"
    if stage == "retirement_income":
        return "从追求收益转向现金流稳定，优先保证提款连续性和回撤可控。"
    return "进一步降低风险资产波动暴露，优先保留流动性和保障性资产。"


def wealth_band_label(band):
    return {
        "low_wealth": "低财富",
        "mid_wealth": "中财富",
        "high_wealth": "高财富",
    }.get(band, band)


def choose_wealth_band(wealth_bands, current_wealth):
    if current_wealth is None or not wealth_bands:
        return "mid_wealth"
    return min(
        wealth_bands.keys(),
        key=lambda name: abs(float(wealth_bands[name]["cash"]) - float(current_wealth)),
    )


def build_advisor_insights(best, fixed, client_profile=None):
    checkpoints = (best.get("lifecycle_checkpoints") or {}).get("checkpoints") or []
    if not checkpoints:
        return None

    client_profile = client_profile or {}
    current_wealth = client_profile.get("current_wealth")
    retirement_age = int(fixed.get("tr", 65))
    lines = []
    for row in checkpoints:
        wealth_bands = row.get("wealth_bands") or {}
        selected_band = choose_wealth_band(wealth_bands, current_wealth)
        selected = wealth_bands.get(selected_band, {
            "alpha": row["mid_wealth_alpha"],
            "consumption": row["mid_wealth_consumption"],
            "cash": None,
        })
        risk_level = describe_risk_level(selected["alpha"])
        stage = phase_label(row["age"], retirement_age)
        portfolio = portfolio_mix_from_alpha(selected["alpha"])
        if row["phase"] == "working":
            client_message = (
                f"客户处于{phase_label_cn(stage)}，按{wealth_band_label(selected_band)}路径，"
                f"建议以{bucket_label_cn(portfolio['bucket'])}为主，当前重点是长期增值与退休前准备。"
            )
        else:
            client_message = (
                f"客户已进入{phase_label_cn(stage)}，按{wealth_band_label(selected_band)}路径，"
                f"建议以{bucket_label_cn(portfolio['bucket'])}承接退休后的现金流需求，减少大幅回撤风险。"
            )
        lines.append({
            "age": row["age"],
            "phase": row["phase"],
            "stage": stage,
            "risk_level": risk_level,
            "portfolio": portfolio,
            "selected_wealth_band": selected_band,
            "selected_band_label": wealth_band_label(selected_band),
            "selected_cash_reference": selected.get("cash"),
            "selected_alpha": selected["alpha"],
            "selected_consumption": selected["consumption"],
            "wealth_bands": wealth_bands,
            "mid_wealth_alpha": row["mid_wealth_alpha"],
            "mid_wealth_consumption": row["mid_wealth_consumption"],
            "portfolio_action": portfolio_action(stage, row["phase"]),
            "client_message": client_message,
        })
    return lines


def render_customer_manager_report(best, report_meta, fixed, client_profile=None):
    checkpoints = (best.get("lifecycle_checkpoints") or {}).get("checkpoints") or []
    if not checkpoints:
        return "# Customer Manager Report\n\n暂无可用于客户沟通的生命周期检查点。\n"

    client_profile = client_profile or {}
    current_wealth = client_profile.get("current_wealth")
    insights = build_advisor_insights(best, fixed, client_profile) or []
    retirement_age = int(fixed.get("tr", 65))
    death_age = int(fixed.get("td", 100))
    objective = best.get("objective", report_meta.get("objective", "maximize_terminal_wealth"))
    params = best.get("params", {})

    first = insights[0]
    last = insights[-1]
    summary = (
        f"客户当前处于`{phase_label_cn(first['stage'])}`，按`{first['selected_band_label']}`路径，"
        f"建议以`{bucket_label_cn(first['portfolio']['bucket'])}`作为主配置思路，"
        f"并在`{retirement_age}岁`退休后逐步转向更稳健的提款型组合。"
    )

    lines = [
        "# Customer Manager Report",
        "",
        "## 客户结论",
        "",
        summary,
        "",
        f"- 生命周期区间: {fixed.get('tb', 20)}岁开始工作，{retirement_age}岁退休，{death_age}岁寿命终点",
        f"- 模型目标: {objective}",
        f"- 最优场景分数: {best.get('score', 0.0):.6f}",
        (
            "- 财富输入: 未提供；本报告会展示低/中/高财富三种路径，并默认按中财富路径生成主结论"
            if current_wealth is None
            else f"- 财富输入: current_wealth={current_wealth}，已按最接近的 policy function 节点取值"
        ),
        f"- 建议起始组合: {bucket_label_cn(first['portfolio']['bucket'])}，权益类约{first['portfolio']['equity_pct']}%，稳健类约{first['portfolio']['stabilizer_pct']}%",
        f"- 模型参数附录: rho={params.get('rho')}, delta={params.get('delta')}, psi={params.get('psi')}, mu={params.get('mu')}, sigr={params.get('sigr')}",
        "",
        "## 推荐组合",
        "",
    ]

    for row in insights:
        lines.append(
            f"- {row['age']}岁 | {phase_label_cn(row['stage'])} | 推荐 `{bucket_label_cn(row['portfolio']['bucket'])}` | "
            f"权益类约{row['portfolio']['equity_pct']}%，稳健类约{row['portfolio']['stabilizer_pct']}% | "
            f"财富路径: {row['selected_band_label']} (参考 cash={row['selected_cash_reference']}) | "
            f"风险级别: {risk_label_cn(row['risk_level'])} | "
            f"原因: {row['portfolio_action']}"
        )

    lines.extend([
        "",
        "## 财富敏感度",
        "",
    ])

    for row in insights:
        wealth_bands = row["wealth_bands"]
        lines.append(
            f"- {row['age']}岁: "
            f"低财富 alpha={wealth_bands['low_wealth']['alpha']:.4f}, "
            f"中财富 alpha={wealth_bands['mid_wealth']['alpha']:.4f}, "
            f"高财富 alpha={wealth_bands['high_wealth']['alpha']:.4f}"
        )

    lines.extend([
        "",
        "## 调仓方向",
        "",
    ])

    for current, nxt in zip(insights, insights[1:]):
        lines.append(
            f"- 从{current['age']}岁到{nxt['age']}岁: 组合从 `{bucket_label_cn(current['portfolio']['bucket'])}` 过渡到 "
            f"`{bucket_label_cn(nxt['portfolio']['bucket'])}`，权益类建议由约{current['portfolio']['equity_pct']}% 调整为约{nxt['portfolio']['equity_pct']}%。"
        )

    lines.extend([
        "",
        "## 客户经理沟通话术",
        "",
    ])

    for row in insights:
        lines.append(f"- {row['age']}岁节点: {row['client_message']}")

    lines.extend([
        "",
        "## 风险提示",
        "",
        f"- 退休年龄是组合切换的关键拐点。如果客户真实退休计划不是 {retirement_age} 岁，建议重新测算。",
        "- 报告给的是配置方向，不是具体基金或产品名单；落地时还要结合客户风险评级和产品准入。",
        "- 如果客户近期有大额支出、医疗安排或提前退休计划，客户经理应人工覆盖模型建议。",
        "",
        "## 下一步动作",
        "",
        "- 先确认客户真实退休年龄、退休后年支出和可接受回撤区间。",
        "- 再确认客户当前 wealth 落在哪条 policy function 路径，不要直接套用中财富结论。",
        "- 最后形成可执行调仓单，明确哪些资产保留、哪些逐步降低、哪些用于退休现金流准备。",
        "",
        "## 关键观察",
        "",
        f"- 退休前主线: {insights[0]['age']}岁到{insights[1]['age']}岁，组合核心仍是 `{bucket_label_cn(insights[0]['portfolio']['bucket'])}`。",
        f"- 退休切换: {retirement_age - 1}岁到{retirement_age}岁，建议开始把沟通重点从增值转到提款稳定。",
        f"- 长寿阶段: {last['age']}岁时，权益类建议降到约{last['portfolio']['equity_pct']}%，更强调保障属性。",
        "",
        f"_生成说明: 共评估 {report_meta.get('total_scenarios')} 个场景，模型运行 {report_meta.get('elapsed_seconds')} 秒。_",
        "",
    ])
    return "\n".join(lines)


def run_real_model(params, fixed, artifact_dir: Path, fast_mode: bool, timeout_sec: int):
    artifact_dir.mkdir(parents=True, exist_ok=True)
    root = repo_root()
    needed = ["f_spline.m", "f_sc_splint.m", "f_ntoil.m", "f_randn.m"]
    for fn in needed:
        shutil.copy2(root / fn, artifact_dir / fn)

    patch_lifecycle_script(root / "life_cycle.m", artifact_dir / "life_cycle.m", params, fixed, fast_mode)

    cmd = ["octave", "--quiet", "life_cycle.m"]
    t0 = time.time()
    proc = subprocess.run(
        cmd,
        cwd=str(artifact_dir),
        capture_output=True,
        text=True,
        encoding="utf-8",
        errors="ignore",
        timeout=timeout_sec,
    )
    dt = time.time() - t0
    (artifact_dir / "octave_stdout.txt").write_text(proc.stdout or "", encoding="utf-8")
    (artifact_dir / "octave_stderr.txt").write_text(proc.stderr or "", encoding="utf-8")
    metric = parse_terminal_metric(artifact_dir)
    year_files = sorted(artifact_dir.glob("year*.txt"))
    policy_summary = build_policy_summary(artifact_dir, fixed) if proc.returncode == 0 and year_files else None
    if policy_summary:
        (artifact_dir / "policy_summary.json").write_text(json.dumps(policy_summary, indent=2, ensure_ascii=False), encoding="utf-8")
        (artifact_dir / "policy_summary.md").write_text(render_policy_markdown(policy_summary), encoding="utf-8")
    lifecycle_checkpoints = build_lifecycle_checkpoints(policy_summary)
    status = "ok" if proc.returncode == 0 and len(year_files) > 0 else "error"
    return {
        "status": status,
        "artifact_dir": str(artifact_dir),
        "code": proc.returncode,
        "metric": metric,
        "run_seconds": round(dt, 3),
        "octave_command": " ".join(cmd),
        "year_files_count": len(year_files),
        "policy_summary_path": str(artifact_dir / "policy_summary.json") if policy_summary else None,
        "lifecycle_checkpoints": lifecycle_checkpoints,
    }


def run_model(params, fixed, artifact_dir: Path, dry_run: bool, use_real_model: bool, fast_mode: bool, timeout_sec: int):
    artifact_dir.mkdir(parents=True, exist_ok=True)
    if dry_run:
        return {"status": "dry_run", "artifact_dir": str(artifact_dir), "metric": None, "run_seconds": 0.0, "year_files_count": 0}
    if use_real_model:
        if not octave_available():
            return {"status": "error_no_octave", "artifact_dir": str(artifact_dir), "metric": None, "run_seconds": 0.0, "year_files_count": 0}
        try:
            return run_real_model(params, fixed, artifact_dir, fast_mode, timeout_sec)
        except subprocess.TimeoutExpired:
            return {"status": "error_timeout", "artifact_dir": str(artifact_dir), "metric": None, "run_seconds": timeout_sec, "year_files_count": 0}
    return {"status": "simulated", "artifact_dir": str(artifact_dir), "metric": None, "run_seconds": 0.0, "year_files_count": 0}


def build_candidates(cfg, method, max_evals, seed):
    keys = ["rho", "delta", "psi", "mu", "sigr"]
    spaces = [cfg["search_space"][k] for k in keys]
    all_points = [dict(zip(keys, tup)) for tup in itertools.product(*spaces)]

    if method == "grid":
        if max_evals and max_evals < len(all_points):
            random.Random(seed).shuffle(all_points)
            return all_points[:max_evals]
        return all_points

    rng = random.Random(seed)
    if not max_evals:
        max_evals = min(200, len(all_points))
    return [rng.choice(all_points) for _ in range(max_evals)]


def parse_args():
    ap = argparse.ArgumentParser(description="LifeCycle optimizer runner")
    ap.add_argument("--config", required=True)
    ap.add_argument("--output-dir", required=True)
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--use-real-model", action="store_true", help="run actual life_cycle.m via octave")
    ap.add_argument("--allow-proxy-fallback", action="store_true", help="when real model fails, keep proxy score instead of exiting")
    ap.add_argument("--fast-mode", action="store_true", help="reduce grid/simulation size for quicker runtime")
    ap.add_argument("--timeout-sec", type=int, default=600)
    ap.add_argument("--search-method", choices=["grid", "random"], default=None)
    ap.add_argument("--max-evals", type=int, default=None)
    ap.add_argument("--progress-every", type=int, default=20)
    return ap.parse_args()


def main():
    args = parse_args()
    cfg = load_config(Path(args.config))
    validate(cfg)

    seed = cfg.get("seed", 0)
    random.seed(seed)
    method = args.search_method or cfg.get("search_method", "grid")
    max_evals = args.max_evals or cfg.get("max_evals")
    fixed = cfg.get("fixed", {})

    out = resolve_output_dir(args.output_dir)
    out.mkdir(parents=True, exist_ok=True)

    objective = cfg.get("objective", "maximize_terminal_wealth")
    candidates = build_candidates(cfg, method, max_evals, seed)
    total = len(candidates)
    if total == 0:
        raise ValueError("no candidates generated")

    if args.use_real_model and not octave_available() and not args.allow_proxy_fallback:
        raise RuntimeError("--use-real-model requires Octave. Install Octave or add --allow-proxy-fallback.")

    octv = shutil.which("octave")
    print(f"[optimizer] method={method} objective={objective} candidates={total} dry_run={args.dry_run} use_real_model={args.use_real_model} fast_mode={args.fast_mode}")
    print(f"[optimizer] octave_path={octv if octv else 'NOT_FOUND'}")
    sys.stdout.flush()

    start = time.time()
    results = []
    results_path = out / "results.jsonl"
    with results_path.open("w", encoding="utf-8") as f:
        for i, params in enumerate(candidates):
            scenario_dir = out / f"scenario_{i:04d}"
            run_meta = run_model(params, fixed, scenario_dir, args.dry_run, args.use_real_model, args.fast_mode, args.timeout_sec)
            if run_meta.get("metric") is not None:
                score = run_meta["metric"]
            elif args.use_real_model and not args.allow_proxy_fallback:
                score = float("-inf")
            else:
                score = score_proxy(params, objective)
            row = {
                "id": i,
                "params": params,
                "objective": objective,
                "score": score,
                "status": run_meta["status"],
                "artifact_dir": run_meta["artifact_dir"],
                "metric": run_meta.get("metric"),
                "run_seconds": run_meta.get("run_seconds", 0.0),
                "year_files_count": run_meta.get("year_files_count", 0),
                "octave_command": run_meta.get("octave_command"),
                "policy_summary_path": run_meta.get("policy_summary_path"),
                "lifecycle_checkpoints": run_meta.get("lifecycle_checkpoints"),
            }
            f.write(json.dumps(row, ensure_ascii=False) + "\n")
            results.append(row)

            if (i + 1) % max(1, args.progress_every) == 0 or i == total - 1:
                elapsed = time.time() - start
                print(f"[optimizer] progress {i+1}/{total} elapsed={elapsed:.1f}s")
                sys.stdout.flush()

    results.sort(key=lambda r: r["score"], reverse=True)
    best = results[0]
    top_k = results[: cfg.get("top_k", 5)]
    (out / "best_params.json").write_text(json.dumps(best, indent=2, ensure_ascii=False), encoding="utf-8")

    failed = sum(1 for r in results if str(r.get("status","")).startswith("error"))
    total_runtime = round(sum(float(r.get("run_seconds", 0.0)) for r in results), 3)

    report = {
        "objective": objective,
        "search_method": method,
        "use_real_model": args.use_real_model,
        "fast_mode": args.fast_mode,
        "total_scenarios": total,
        "failed_scenarios": failed,
        "total_model_runtime_seconds": total_runtime,
        "elapsed_seconds": round(time.time() - start, 3),
        "best": best,
        "top_k": top_k,
    }
    (out / "report.md").write_text("# LifeCycle Optimization Report\n\n" + json.dumps(report, indent=2, ensure_ascii=False), encoding="utf-8")
    (out / "customer_manager_report.md").write_text(
        render_customer_manager_report(best, report, fixed, cfg.get("client_profile")),
        encoding="utf-8",
    )
    print(f"[optimizer] done best_score={best['score']:.6f} output_dir={out}")
    if args.use_real_model and failed > 0:
        print(f"[optimizer] warning: {failed} real-model scenarios failed; check scenario logs")
    print(f"[optimizer] absolute_output_dir={out.resolve()}")


def optimize(config_path, output_dir,
             use_real_model=False,
             fast_mode=False,
             dry_run=False):
    """
    Callable version of optimizer (for skill / agent use)
    """

    # ===== 模拟 argparse =====
    class Args:
        pass

    args = Args()
    args.config = config_path
    args.output_dir = output_dir
    args.dry_run = dry_run
    args.use_real_model = use_real_model
    args.allow_proxy_fallback = True
    args.fast_mode = fast_mode
    args.timeout_sec = 600
    args.search_method = None
    args.max_evals = None
    args.progress_every = 20

    # ===== 下面就是 main() 的逻辑（轻微改造） =====
    cfg = load_config(Path(args.config))
    validate(cfg)

    seed = cfg.get("seed", 0)
    random.seed(seed)
    method = args.search_method or cfg.get("search_method", "grid")
    max_evals = args.max_evals or cfg.get("max_evals")
    fixed = cfg.get("fixed", {})

    out = resolve_output_dir(args.output_dir)
    out.mkdir(parents=True, exist_ok=True)

    objective = cfg.get("objective", "maximize_terminal_wealth")
    candidates = build_candidates(cfg, method, max_evals, seed)
    total = len(candidates)

    start = time.time()
    results = []

    for i, params in enumerate(candidates):
        scenario_dir = out / f"scenario_{i:04d}"

        run_meta = run_model(
            params, fixed, scenario_dir,
            args.dry_run, args.use_real_model,
            args.fast_mode, args.timeout_sec
        )

        if run_meta.get("metric") is not None:
            score = run_meta["metric"]
        else:
            score = score_proxy(params, objective)

        row = {
            "id": i,
            "params": params,
            "score": score,
            "status": run_meta["status"],
            "lifecycle_checkpoints": run_meta.get("lifecycle_checkpoints"),
        }
        results.append(row)

    # ===== 排序 + best =====
    results.sort(key=lambda r: r["score"], reverse=True)
    best = results[0]

    return {
        "best": best,
        "top_k": results[:5],
        "total": total,
        "runtime": round(time.time() - start, 3),
        "output_dir": str(out)
    }


def _safe_load_policy_summary(best):
    path = best.get("policy_summary_path")
    if not path:
        return None
    p = Path(path)
    if not p.exists():
        return None
    with p.open("r", encoding="utf-8") as f:
        return json.load(f)


def _band_point(row, band):
    wealth_bands = row.get("wealth_bands")
    if wealth_bands and band in wealth_bands:
        return wealth_bands[band]
    return {
        "alpha": row.get("mid_wealth_alpha"),
        "consumption": row.get("mid_wealth_consumption"),
        "cash": None,
    }


def _build_decade_rows(best, fixed, client_profile):
    policy_summary = _safe_load_policy_summary(best)
    current_wealth = (client_profile or {}).get("current_wealth")
    current_age = (client_profile or {}).get("current_age", fixed.get("tb", 20))
    rows = []

    if policy_summary and policy_summary.get("per_age"):
        per_age = policy_summary["per_age"]
        age_map = {row["age"]: row for row in per_age}
        ages = list(range(int(fixed.get("tb", 20)), int(fixed.get("td", 100)), 10))
        for start_age in ages:
            src = age_map.get(start_age)
            if not src:
                nearest = next((row for row in per_age if row["age"] >= start_age), None)
                if not nearest:
                    continue
                src = nearest
            wealth_bands = {
                band: {
                    "alpha": src["alpha"]["points"][band]["value"],
                    "consumption": src["consumption"]["points"][band]["value"],
                    "cash": src["alpha"]["points"][band]["cash"],
                }
                for band in ["low_wealth", "mid_wealth", "high_wealth"]
            }
            selected_band = choose_wealth_band(wealth_bands, current_wealth)
            selected = wealth_bands[selected_band]
            portfolio = portfolio_mix_from_alpha(selected["alpha"])
            rows.append({
                "start_age": start_age,
                "end_age": min(start_age + 10, int(fixed.get("td", 100))),
                "selected_band": selected_band,
                "selected_band_label": wealth_band_label(selected_band),
                "selected_cash": selected["cash"],
                "alpha": selected["alpha"],
                "consumption": selected["consumption"],
                "portfolio": portfolio,
                "wealth_bands": wealth_bands,
                "is_current_window": start_age <= current_age < min(start_age + 10, int(fixed.get("td", 100))),
            })
        return rows

    checkpoints = (best.get("lifecycle_checkpoints") or {}).get("checkpoints") or []
    for i, row in enumerate(checkpoints):
        wealth_bands = row.get("wealth_bands") or {
            "mid_wealth": {
                "alpha": row.get("mid_wealth_alpha"),
                "consumption": row.get("mid_wealth_consumption"),
                "cash": None,
            }
        }
        selected_band = choose_wealth_band(wealth_bands, current_wealth)
        selected = _band_point(row, selected_band)
        portfolio = portfolio_mix_from_alpha(selected["alpha"])
        start_age = row["age"]
        end_age = checkpoints[i + 1]["age"] if i + 1 < len(checkpoints) else int(fixed.get("td", 100))
        rows.append({
            "start_age": start_age,
            "end_age": end_age,
            "selected_band": selected_band,
            "selected_band_label": wealth_band_label(selected_band),
            "selected_cash": selected["cash"],
            "alpha": selected["alpha"],
            "consumption": selected["consumption"],
            "portfolio": portfolio,
            "wealth_bands": wealth_bands,
            "is_current_window": start_age <= current_age < end_age,
        })
    return rows


def render_customer_manager_report(best, report_meta, fixed, client_profile=None):
    client_profile = client_profile or {}
    current_age = client_profile.get("current_age", fixed.get("tb", 20))
    current_wealth = client_profile.get("current_wealth")
    rows = _build_decade_rows(best, fixed, client_profile)
    if not rows:
        return "# Customer Manager Report\n\n暂无可用于客户沟通的建议。\n"

    current_row = next((row for row in rows if row["is_current_window"]), rows[0])
    params = best.get("params", {})
    lines = [
        "# Customer Manager Report",
        "",
        "## 客户当前信息",
        "",
        f"- 当前年龄: {current_age}岁",
        f"- 当前财富: {'未提供' if current_wealth is None else f'current_wealth={current_wealth}'}",
        f"- 当前命中财富路径: {current_row['selected_band_label']}",
        f"- 当前建议组合: {bucket_label_cn(current_row['portfolio']['bucket'])}",
        f"- 当前建议配比: 权益类约{current_row['portfolio']['equity_pct']}%，稳健类约{current_row['portfolio']['stabilizer_pct']}%",
        f"- 当前建议消费强度参考: {current_row['consumption']:.4f}",
        f"- 模型目标: {best.get('objective', report_meta.get('objective', 'maximize_terminal_wealth'))}",
        f"- 最优场景分数: {best.get('score', 0.0):.6f}",
        f"- 参数附录: rho={params.get('rho')}, delta={params.get('delta')}, psi={params.get('psi')}, mu={params.get('mu')}, sigr={params.get('sigr')}",
        "",
        "## 后续建议",
        "",
    ]

    for row in rows:
        wealth_bands = row["wealth_bands"]
        lines.append(
            f"- {row['start_age']}岁到{row['end_age']}岁: 建议以`{bucket_label_cn(row['portfolio']['bucket'])}`为主，"
            f"权益类约{row['portfolio']['equity_pct']}%，稳健类约{row['portfolio']['stabilizer_pct']}%。"
        )
        lines.append(
            f"  财富路径: {row['selected_band_label']} (参考 cash={row['selected_cash']})；"
            f"低/中/高财富 alpha 分别为 "
            f"{wealth_bands.get('low_wealth', {}).get('alpha', row['alpha']):.4f}/"
            f"{wealth_bands.get('mid_wealth', {}).get('alpha', row['alpha']):.4f}/"
            f"{wealth_bands.get('high_wealth', {}).get('alpha', row['alpha']):.4f}。"
        )
        lines.append(
            f"  建议说明: {portfolio_action('mid_accumulation', 'working') if row['start_age'] < fixed.get('tr', 65) else '继续降低波动，强化流动性和稳定性。'}"
        )

    lines.extend([
        "",
        "## 客户经理使用提示",
        "",
        "- 先确认客户当前年龄和财富水平，再选择对应的 decade 建议，不要直接套用中财富路径。",
        "- 如果客户财富明显高于或低于当前路径参考 cash，优先改用更接近的财富档位。",
        "- 报告给的是配置方向，不是具体产品清单；落地时还要结合客户风险评级和产品准入。",
        "",
        f"_生成说明: 共评估 {report_meta.get('total_scenarios')} 个场景，模型运行 {report_meta.get('elapsed_seconds')} 秒。_",
        "",
    ])
    return "\n".join(lines)


if __name__ == "__main__":
    main()
