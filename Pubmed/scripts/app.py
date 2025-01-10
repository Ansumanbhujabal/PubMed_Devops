from fastapi import FastAPI
from fastapi.middleware.cors import CORSMiddleware
from Bio import Entrez
import pandas as pd
import re
from typing import List

# Initialize FastAPI app
app = FastAPI()

# Enable CORS for Streamlit to access the API
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_methods=["*"],
    allow_headers=["*"],
)

# Set email for NCBI
Entrez.email = "your_email@example.com"

# Helper functions for PubMed search and processing
def search(query: str) -> List[str]:
    handle = Entrez.esearch(db="pubmed", sort="relevance", retmax="30", retmode="xml", term=query)
    results = Entrez.read(handle)
    return results["IdList"]

def fetch_details(id_list: List[str]) -> List:
    ids = ",".join(id_list)
    handle = Entrez.efetch(db="pubmed", retmode="xml", id=ids)
    return Entrez.read(handle)["PubmedArticle"]

def extract_email(affiliation: str) -> str:
    match = re.search(r"[\w\.-]+@[\w\.-]+", affiliation)
    return match.group(0) if match else ""

def identify_academic_non_academic_authors(authors: List[str], affiliations: List[str]) -> List[str]:
    non_academic_authors = []
    academic_keywords = ["school", "college", "university", "institute", "department", "research"]
    for author, affiliation in zip(authors, affiliations):
        if not any(keyword in affiliation.lower() for keyword in academic_keywords):
            non_academic_authors.append(f"{author} ({affiliation})")
    return non_academic_authors

def process_studies(query: str, company_keywords: List[str]) -> pd.DataFrame:
    id_list = search(query)
    papers = fetch_details(id_list)
    data = []


    for paper in papers:
        paper_id = paper["MedlineCitation"]["PMID"]
        title = paper["MedlineCitation"]["Article"]["ArticleTitle"]
        try:
            abstract = paper["MedlineCitation"]["Article"]["Abstract"]["AbstractText"][0]
        except KeyError:
            abstract = "No Abstract"
        journal = paper["MedlineCitation"]["Article"]["Journal"]["Title"]
        language = paper["MedlineCitation"]["Article"]["Language"][0]
        authors, affiliations,all_affiliations = [], [],[]
        for author in paper["MedlineCitation"]["Article"]["AuthorList"]:
            try:
                authors.append(author["LastName"] + " " + author["ForeName"])
                if "AffiliationInfo" in author and author["AffiliationInfo"]:
                    affiliations.append(author["AffiliationInfo"][0]["Affiliation"])
            except KeyError:
                pass
        all_affiliations = " ".join(affiliations)
        email = extract_email(all_affiliations)                
        non_academic_authors = identify_academic_non_academic_authors(authors, affiliations)
        if non_academic_authors:
            data.append({
                "ID": paper_id,
                "Title": title,
                "Abstract": abstract,
                "Journal": journal,
                "Language": language,
                "Name":"; ".join(authors),
                "Mail":email,
                "Non-academic Authors": "; ".join(non_academic_authors)
            })
    return pd.DataFrame(data)

@app.post("/process/")
def process(query: str, keywords: str, operator: str):
    # Modify query based on operator and keywords
    if keywords:
        keyword_list = keywords.split(",")
        keyword_query = f" {operator} ".join(keyword_list)
        query = f"({query}) {operator} ({keyword_query})"
    df = process_studies(query, keywords.split(","))
    return df.to_dict(orient="records")
