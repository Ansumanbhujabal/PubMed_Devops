import streamlit as st
import pandas as pd
import requests

# FastAPI backend URL
BACKEND_URL = "http://127.0.0.1:8080/process/"

# Streamlit UI
st.title("PubMed Search and Analysis")

# Input fields
query = st.text_input("Query", value="CHILDREN AND EPILEPSY")
keywords = st.text_input("Keywords (comma-separated)", value="pharmaceutical,biotech,bioscience")
operator = st.selectbox("Logical Operator", ["AND", "OR", "NOT"])

# Process button
if st.button("Process"):
    params = {"query": query, "keywords": keywords, "operator": operator}
    response = requests.post(BACKEND_URL, params=params)
    if response.status_code == 200:
        data = response.json()
        if data:
            st.write("Filtered Results:")
            df = pd.DataFrame(data)
            st.dataframe(df)
            csv = df.to_csv(index=False).encode("utf-8")
            st.download_button(
                label="Download CSV",
                data=csv,
                file_name="filtered_results.csv",
                mime="text/csv"
            )
        else:
            st.write("No results found.")
    else:
        st.write("Error fetching results.")

