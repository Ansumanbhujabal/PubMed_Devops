import requests
import csv
import re
from bs4 import BeautifulSoup
from typing import List, Dict

# PubMed API Base URL
PUBMED_API_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

# Companies Keywords List (can be expanded)
COMPANY_KEYWORDS = [
    "pharmaceutical", "biotech", "bioscience", "life sciences", "drug discovery"
]

def fetch_pubmed_ids(query: str) -> List[str]:
    """
    Fetch PubMed IDs based on the query.
    """
    search_url = f"{PUBMED_API_BASE}/esearch.fcgi"
    params = {
        "db": "pubmed",
        "term": query,
        "retmax": 100,  # Adjust based on limits
        "retmode": "json",
    }
    response = requests.get(search_url, params=params)
    response.raise_for_status()
    data = response.json()
    return data.get("esearchresult", {}).get("idlist", [])

def fetch_paper_details(pubmed_id: str) -> Dict:
    """
    Fetch detailed information about a PubMed paper.
    """
    fetch_url = f"{PUBMED_API_BASE}/efetch.fcgi"
    params = {
        "db": "pubmed",
        "id": pubmed_id,
        "retmode": "xml",
    }
    response = requests.get(fetch_url, params=params)
    response.raise_for_status()
    soup = BeautifulSoup(response.content, "xml")

    title = soup.find("ArticleTitle").text if soup.find("ArticleTitle") else ""
    pub_date = soup.find("PubDate").text if soup.find("PubDate") else ""
    authors = [author.find("LastName").text + " " + author.find("ForeName").text
               for author in soup.find_all("Author")]
    affiliations = [aff.text for aff in soup.find_all("Affiliation")]

    # Extract corresponding author email and company affiliations
    corresponding_email = extract_email(" ".join(affiliations))
    company_affiliations, non_academic_authors = extract_company_affiliations(authors, affiliations)

    return {
        "PubmedID": pubmed_id,
        "Title": title,
        "Publication Date": pub_date,
        "Non-academic Author(s)": ", ".join(non_academic_authors),
        "Company Affiliation(s)": ", ".join(company_affiliations),
        "Corresponding Author Email": corresponding_email,
    }

def extract_email(affiliation_text: str) -> str:
    """
    Extract email address from the affiliation text.
    """
    email_match = re.search(r"[\w\.-]+@[\w\.-]+", affiliation_text)
    return email_match.group(0) if email_match else ""

def extract_company_affiliations(authors: List[str], affiliations: List[str]) -> (List[str], List[str]):
    """
    Identify authors with company affiliations.
    """
    company_affiliations = []
    non_academic_authors = []
    for author, affiliation in zip(authors, affiliations):
        if any(keyword in affiliation.lower() for keyword in COMPANY_KEYWORDS):
            company_affiliations.append(affiliation)
            non_academic_authors.append(author)
    return company_affiliations, non_academic_authors

def save_to_csv(data: List[Dict], filename: str) -> None:
    """
    Save the fetched data to a CSV file.
    """
    with open(filename, mode="w", newline="", encoding="utf-8") as file:
        writer = csv.DictWriter(file, fieldnames=data[0].keys())
        writer.writeheader()
        writer.writerows(data)

def main(query: str, output_file: str):
    """
    Main function to fetch and process papers.
    """
    pubmed_ids = fetch_pubmed_ids(query)
    results = []
    for pubmed_id in pubmed_ids:
        details = fetch_paper_details(pubmed_id)
        results.append(details)

    save_to_csv(results, output_file)
    print(f"Results saved to {output_file}")

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Fetch research papers from PubMed.")
    parser.add_argument("query", type=str, help="PubMed query string.")
    parser.add_argument("output", type=str, help="Output CSV file name.")

    args = parser.parse_args()
    main(args.query, args.output)
