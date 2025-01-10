
---

# PubMed Search GUI Module

This project provides an API for searching and analyzing PubMed articles, using FastAPI for the backend and Streamlit for the frontend. It fetches PubMed data based on user queries, processes the results, and allows users to download the filtered data in CSV format.

## Features
- **FastAPI Backend**: Fetches PubMed articles, processes data (e.g., extracting author details, affiliations), and allows filtering based on keywords and logical operators.
- **Streamlit Frontend**: Provides a simple UI to input queries, define keywords, and display results in a tabular format.
- **Email Extraction**: Automatically extracts email addresses from affiliations.
- **Non-Academic Author Identification**: Identifies and displays non-academic authors based on their affiliations.
- **CSV Export**: Allows users to download filtered results as a CSV file.

## Technologies Used
- **FastAPI**: For building the backend API to interact with PubMed.
- **Streamlit**: For creating the frontend UI.
- **BioPython**: For interacting with the NCBI Entrez API.
- **pandas**: For managing and processing the data.

## Setup Instructions

### 1. Clone the repository
```bash
git clone https://github.com/Ansumanbhujabal/PubMed_Devops.git
cd PubMed
```

### 2. Install dependencies
This project uses Poetry for dependency management.

```bash
poetry install
```

### 3. Set up your NCBI Email
To interact with PubMed, you need to set up your email in the FastAPI backend. Edit the following line in `app.py`:

```python
Entrez.email = "your_email@example.com"
```

### 4. Run the Application
To start both the FastAPI backend and the Streamlit frontend:

```bash
get-papers-list
```

This will start the FastAPI backend at `http://127.0.0.1:8080/` and the Streamlit frontend will be available at the default port (usually `8501`).

### 5. Usage
1. **Backend**: You can access the API at `http://127.0.0.1:8080/process/` for querying PubMed articles.
2. **Frontend**: Go to the Streamlit UI in your browser and enter your query, keywords, and logical operator to filter the results.

### 6. Exporting Results
Once you get the filtered results, you can download them as a CSV by clicking the "Download CSV" button on the Streamlit UI.

## Example Query

**Query**: `CHILDREN AND EPILEPSY`  
**Keywords**: `pharmaceutical,biotech,bioscience`  
**Logical Operator**: `AND`

This query will return PubMed articles related to "Children and Epilepsy" and filter them by the provided keywords.

## Project Structure
- `scripts/`: Contains the main scripts for running the FastAPI backend and Streamlit frontend.
- `ui.py`: Streamlit frontend code.
- `app.py`: FastAPI backend code.
- `pyproject.toml`: The Poetry configuration file for managing dependencies.

## Dependencies
The project requires the following Python packages:

- `argparse >=1.4.0,<2.0.0`
- `logging >=0.4.9.6,<0.5.0.0`
- `bio >=1.7.1,<2.0.0`
- `pandas >=2.2.3,<3.0.0`
- `tabulate >=0.9.0,<0.10.0`
- `colorama >=0.4.6,<0.5.0`
- `streamlit >=1.41.1,<2.0.0`
- `requests >=2.32.3,<3.0.0`

## License
This project is licensed under the MIT License.
## Docs
https://docs.google.com/document/d/1ieJvURbmxlTRyXHsCCEMCUGauku0kFwniOCrXv3daFw/edit?usp=sharing
https://docs.google.com/document/d/1QDn6nmw46fr7FCWWJNJOn3UEE0TB9xgVL0vtzSKZdfs/edit?usp=sharing
## Author
- **Ansuman SS Bhujabala**

---
