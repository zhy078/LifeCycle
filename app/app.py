import sys
from pathlib import Path

import streamlit as st


APP_DIR = Path(__file__).resolve().parent
if str(APP_DIR) not in sys.path:
    sys.path.insert(0, str(APP_DIR))

from advisor import DEFAULT_RESULTS_DIR, build_advice


st.set_page_config(page_title="Lifecycle Advisor", layout="wide")
st.title("Lifecycle Advisor")

client_text = st.text_input("Client input", value="45岁，100万，偏稳健")
results_dir = st.text_input("Completed results directory", value=str(DEFAULT_RESULTS_DIR))

if st.button("Generate advice", type="primary"):
    result = build_advice(client_text, results_dir=results_dir)
    if result["ok"]:
        st.success("Advice generated")
    else:
        st.error(result["error"])
    st.markdown(result["markdown"])
