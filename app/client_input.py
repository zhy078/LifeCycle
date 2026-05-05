import re


RISK_MAP = {
    "conservative": 10,
    "balanced": 8,
    "aggressive": 6,
}


RISK_KEYWORDS = [
    ("conservative", ["conservative", "stable", "safe", "low risk", "稳健", "保守", "稳定", "低风险", "偏稳健"]),
    ("aggressive", ["aggressive", "growth", "high return", "risky", "激进", "成长", "高收益", "愿意冒险"]),
    ("balanced", ["balanced", "moderate", "neutral", "平衡", "均衡", "中等", "适中"]),
]


def map_risk_preference(text):
    lowered = text.lower()
    for label, keywords in RISK_KEYWORDS:
        if any(keyword in lowered for keyword in keywords):
            return label, RISK_MAP[label]
    return "balanced", RISK_MAP["balanced"]


def _parse_age(text):
    patterns = [
        r"(\d{1,3})\s*岁",
        r"age\s*[:=]?\s*(\d{1,3})",
        r"(\d{1,3})\s*years?\s*old",
    ]
    for pattern in patterns:
        match = re.search(pattern, text, flags=re.IGNORECASE)
        if match:
            return int(match.group(1))
    return None


def _parse_wealth(text):
    text = text.replace(",", "")
    wan = re.search(r"(\d+(?:\.\d+)?)\s*万", text)
    if wan:
        # Model wealth unit convention used in demos: 1 unit = 100k currency.
        return float(wan.group(1)) / 10.0

    million = re.search(r"(\d+(?:\.\d+)?)\s*(million|m)\b", text, flags=re.IGNORECASE)
    if million:
        # Assume 1 model unit = 100k, so 1 million = 10 model units.
        return float(million.group(1)) * 10.0

    thousand = re.search(r"(\d+(?:\.\d+)?)\s*(thousand|thousands|k)\b", text, flags=re.IGNORECASE)
    if thousand:
        # Assume 1 model unit = 100k, so 100 thousand = 1 model unit.
        return float(thousand.group(1)) / 100.0

    wealth = re.search(r"(wealth|资产|财富|本金)\s*[:=]?\s*(\d+(?:\.\d+)?)", text, flags=re.IGNORECASE)
    if wealth:
        return float(wealth.group(2))
    return None


def _parse_income_profile(text):
    lowered = text.lower()
    if any(token in lowered for token in ["retired", "retirement", "no income", "pension only"]):
        return "retired"
    if any(token in lowered for token in ["volatile income", "unstable income", "business income", "freelance", "self-employed"]):
        return "volatile"
    if any(token in lowered for token in ["stable income", "salary", "salaried", "employed", "wage"]):
        return "stable"
    return "stable"


def parse_client_input(text):
    risk_label, rho = map_risk_preference(text)
    return {
        "age": _parse_age(text),
        "wealth": _parse_wealth(text),
        "income_profile": _parse_income_profile(text),
        "risk_preference": risk_label,
        "rho": rho,
        "raw_text": text,
    }
