import argparse
import logging
import re
from Bio import Entrez
import pandas as pd
from tabulate import tabulate
from typing import List, Tuple
from colorama import Fore, Style, init

# Initialize colorama
init(autoreset=True)

# Set email for NCBI
Entrez.email = "your_email@example.com"

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
    handlers=[
        logging.StreamHandler(),
        logging.FileHandler("pubmed_processing.log")
    ]
)
logger = logging.getLogger(__name__)

# Search function to fetch study IDs
def search(query: str) -> List[str]:
    handle = Entrez.esearch(db="pubmed", sort="relevance", retmax="300", retmode="xml", term=query)
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

# Identify academic and non-academic authors
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
    year_list, month_list, date_list = [], [], []
    author_names_list, email_list, all_affiliations_list = [], [], []
    academic_authors_list, non_academic_authors_list = [], []

    for chunk_i in range(0, len(id_list), chunk_size):
        chunk = id_list[chunk_i:chunk_i + chunk_size]
        papers = fetch_details(chunk)

        for paper in papers:
            # Extract ID
            paper_id = paper["MedlineCitation"]["PMID"]
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
            for author in paper["MedlineCitation"]["Article"]["AuthorList"]:
                try:
                    author_name = author["LastName"] + " " + author["ForeName"]
                    authors.append(author_name)
                    if "AffiliationInfo" in author and author["AffiliationInfo"]:
                        affiliations.append(author["AffiliationInfo"][0]["Affiliation"])
                except KeyError:
                    pass

            # Identify corresponding author email
            all_affiliations = " ".join(affiliations)
            email = extract_email(all_affiliations)

            # Identify academic and non-academic authors
            academic_authors, non_academic_authors = identify_academic_non_academic_authors(authors, affiliations)

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
        "Date": date_list,
        "Author Names": author_names_list,
        "Corresponding Author Email": email_list,
        "All Affiliations": all_affiliations_list,
        "Academic Authors": academic_authors_list,
        "Non-academic Authors": non_academic_authors_list,
    })

# Log rows with non-academic authors and save to a filtered CSV
def log_and_save_non_academic_authors(df: pd.DataFrame, output_file: str):
    # Filter rows with non-academic authors
    filtered_df = df[df["Non-academic Authors"] != ""]
    
    if not filtered_df.empty:
        # Log the filtered rows
        summary = filtered_df[["ID", "Journal", "Language", "Year", "Month", "Date", "Author Names", "Corresponding Author Email", "Non-academic Authors"]]
        logger.info(Fore.YELLOW + "\n" + tabulate(summary, headers="keys", tablefmt="grid"))
        
        # Save filtered rows to CSV
        filtered_df.to_csv(output_file, index=False)
        logger.info(Fore.GREEN + f"Filtered results saved to {output_file}")
    else:
        logger.info(Fore.RED + "No rows with non-academic authors found.")

# Main function
def main():
    header_text = """    

         _      _                  _              _   _         _          _                   _             _              _     
        /\ \   /\_\               / /\           /\_\/\_\ _    /\ \       /\ \               /\ \           _\ \           /\ \   
       /  \ \ / / /         _    / /  \         / / / / //\_\ /  \ \     /  \ \____         /  \ \         /\__ \          \ \ \  
      / /\ \ \\ \ \__      /\_\ / / /\ \       /\ \/ \ \/ / // /\ \ \   / /\ \_____\       / /\ \ \       / /_ \_\         /\ \_\ 
     / / /\ \_\\ \___\    / / // / /\ \ \     /  \____\__/ // / /\ \_\ / / /\/___  /      / / /\ \ \     / / /\/_/        / /\/_/ 
    / / /_/ / / \__  /   / / // / /\ \_\ \   / /\/________// /_/_ \/_// / /   / / /      / / /  \ \_\   / / /            / / /    
   / / /__\/ /  / / /   / / // / /\ \ \___\ / / /\/_// / // /____/\  / / /   / / /      / / /    \/_/  / / /            / / /     
  / / /_____/  / / /   / / // / /  \ \ \__// / /    / / // /\____\/ / / /   / / /      / / /          / / / ____       / / /      
 / / /        / / /___/ / // / /____\_\ \ / / /    / / // / /______ \ \ \__/ / /      / / /________  / /_/_/ ___/\ ___/ / /__     
/ / /        / / /____\/ // / /__________\\/_/    / / // / /_______\ \ \___\/ /      / / /_________\/_______/\__\//\__\/_/___\    
\/_/         \/_________/ \/_____________/        \/_/ \/__________/  \/_____/       \/____________/\_______\/    \/_________/    
                                                                                                                                  

M̳a̳d̳e̳ B̳y̳ @̳A̳n̳s̳u̳m̳a̳n̳B̳h̳u̳j̳a̳b̳a̳l̳a̳
""" 

        # Argument parser
    parser = argparse.ArgumentParser(description="PubMed Search CLI")
    parser.add_argument("-q", "--query", type=str, required=True, help="Search query for PubMed")
    parser.add_argument("-o", "--output", type=str, required=True, help="Output CSV file name")
    args = parser.parse_args()

    # Print greeting message
    print(Fore.GREEN + header_text)

    # Fetch query and output file
    query = args.query
    output_file = args.output
    company_keywords = ["pharmaceutical", "biotech", "bioscience", "life sciences", "drug discovery"]

    print(Fore.GREEN + "Searching PubMed...")
    id_list = search(query)

    print(Fore.GREEN + f"Found {len(id_list)} studies. Processing...")
    df = process_studies(id_list, company_keywords)

    print(Fore.GREEN + "Logging and saving rows with non-academic authors...")
    log_and_save_non_academic_authors(df, output_file)

    print(Fore.GREEN + f"Results saved to {output_file}")

if __name__ == "__main__":
    main()

















