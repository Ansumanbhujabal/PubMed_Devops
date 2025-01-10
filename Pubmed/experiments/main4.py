from Bio import Entrez
import pandas as pd
import re
from typing import List, Tuple

# Set email for NCBI
Entrez.email = "your_email@example.com"

# Search function to fetch study IDs
def search(query: str) -> List[str]:
    handle = Entrez.esearch(db="pubmed", sort="relevance", retmax="30", retmode="xml", term=query)
    results = Entrez.read(handle)
    return results["IdList"]

# Fetch details of studies
def fetch_details(id_list: List[str]) -> List:
    ids = ",".join(id_list)
    handle = Entrez.efetch(db="pubmed", retmode="xml", id=ids)
    return Entrez.read(handle)["PubmedArticle"]

# Helper function to extract email from affiliation
def extract_email(affiliation: str) -> str:
    match = re.search(r"[\w\.-]+@[\w\.-]+", affiliation)
    return match.group(0) if match else ""

def identify_academic_non_academic_authors(authors: List[str], affiliations: List[str]) -> Tuple[List[str], List[str]]:
    academic_authors = []
    non_academic_authors = []

    # List of keywords to identify academic affiliations
    academic_keywords = ["school", "college", "university", "institute", "department", "research"]

    for author, affiliation in zip(authors, affiliations):
        if any(keyword in affiliation.lower() for keyword in academic_keywords):
            academic_authors.append(f"{author} ({affiliation})")
        else:
            non_academic_authors.append(f"{author} ({affiliation})")

    return academic_authors, non_academic_authors

# Process studies to create a DataFrame
def process_studies(id_list: List[str], company_keywords: List[str], chunk_size: int = 10000) -> pd.DataFrame:
    id_col, title_list, abstract_list, journal_list, language_list = [], [], [], [], []
    year_list, month_list ,date_list= [], [],[]
    author_names_list, email_list, all_affiliations_list,academic_authors_list, non_academic_authors_list = [], [], [], [],[]

    for chunk_i in range(0, len(id_list), chunk_size):
        chunk = id_list[chunk_i:chunk_i + chunk_size]
        papers = fetch_details(chunk)

        for paper in papers:
            # Extract ID
            paper_id = paper["MedlineCitation"]["PMID"]
            print("------------------------------->>>>>>>>>>>>>>")
            print(f"Paper id is {paper_id}")
            print("------------------------------->>>>>>>>>>>>>>")
            id_col.append(paper_id)
    
            # Extracting publication metadata
            title = paper["MedlineCitation"]["Article"]["ArticleTitle"]
            try:
                abstract = paper["MedlineCitation"]["Article"]["Abstract"]["AbstractText"][0]
            except KeyError:
                abstract = "No Abstract"
            journal = paper["MedlineCitation"]["Article"]["Journal"]["Title"]
            language = paper["MedlineCitation"]["Article"]["Language"][0]
            try:
                year = paper["MedlineCitation"]["Article"]["Journal"]["JournalIssue"]["PubDate"]["Year"]
            except KeyError:
                year = "No Data"
            try:
                month = paper["MedlineCitation"]["Article"]["Journal"]["JournalIssue"]["PubDate"]["Month"]
            except KeyError:
                month = "No Data"
            try:
                date = paper["MedlineCitation"]["Article"]["Journal"]["JournalIssue"]["PubDate"]["Day"]
            except KeyError:
                date = "No Data"                

            # Extracting author and affiliation information
            authors = []
            affiliations = []
            emails = []
            for author in paper["MedlineCitation"]["Article"]["AuthorList"]:
                try:
                    author_name = author["LastName"] + " " + author["ForeName"]
                    authors.append(author_name)
                except KeyError:
                    pass
            if "AffiliationInfo" in author and author["AffiliationInfo"]:
                if "Affiliation" in author["AffiliationInfo"][0]:
                    affiliations.append(author["AffiliationInfo"][0]["Affiliation"])
                  # Extract email if available
                if "Email" in author["AffiliationInfo"][0]:
                    emails.append(author["AffiliationInfo"][0]["Email"])   

            # Identify corresponding author email
            all_affiliations = " ".join(affiliations)
            email = extract_email(all_affiliations)

            # Append data to respective lists
            title_list.append(title)
            abstract_list.append(abstract)
            journal_list.append(journal)
            language_list.append(language)
            year_list.append(year)
            month_list.append(month)
            date_list.append(date)
            author_names_list.append("; ".join(authors))
            email_list.append(email)
            all_affiliations_list.append(all_affiliations)
            academic_authors, non_academic_authors = identify_academic_non_academic_authors( authors,all_affiliations_list)
            academic_authors_list.append("; ".join(academic_authors))
            non_academic_authors_list.append("; ".join(non_academic_authors))            

    # Create DataFrame
    return pd.DataFrame({
        "ID": id_col,
        "Title": title_list,
        "Abstract": abstract_list,
        "Journal": journal_list,
        "Language": language_list,
        "Year": year_list,
        "Month": month_list,
        "Date":date_list,
        "Author Names": author_names_list,
        "Corresponding Author Email": email_list,
        "All Affliations":all_affiliations_list,
        "Academic Authors": academic_authors_list,
        "Non-academic Authors": non_academic_authors_list,
    })


    # Sort by presence of non-academic authors
    df["Has Non-academic Authors"] = df["Non-academic Authors"].apply(lambda x: 1 if x else 0)
    df.sort_values(by="Has Non-academic Authors", ascending=False, inplace=True)
    df.drop(columns="Has Non-academic Authors", inplace=True)

    return df

# Main function
def main():
    query = "CHILDREN AND EPILEPSY"  # Example query
    output_file = "pubmed_results.csv"
    company_keywords = ["pharmaceutical", "biotech", "bioscience", "life sciences", "drug discovery"]

    print("Searching PubMed...")
    id_list = search(query)

    print(f"Found {len(id_list)} studies. Processing...")
    df = process_studies(id_list, company_keywords)

    print("Saving to CSV...")
    df.to_csv(output_file, index=False)
    print(f"Results saved to {output_file}")

if __name__ == "__main__":
    main()
