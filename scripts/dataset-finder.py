#!/usr/bin/env python
# import
## batteries
import os
import sys
import re
import time
import json
import argparse
from pprint import pprint
from datetime import datetime, timedelta
from typing import List, Dict, Literal, Any
import xml.etree.ElementTree as ET
## 3rd party
from Bio import Entrez

# globals
DatabaseType = Literal["sra", "gds"]
start_date = datetime.now() - timedelta(days=7)

# argparse
class CustomFormatter(argparse.ArgumentDefaultsHelpFormatter,
                      argparse.RawDescriptionHelpFormatter):
    pass

desc = 'Find SRA and GEO datasets from NCBI.'
epi = """DESCRIPTION:
This script searches the NCBI Sequence Read Archive (SRA) and Gene Expression Omnibus (GEO)
databases for datasets matching the specified criteria. The script will output JSON files
containing the dataset information.
"""
parser = argparse.ArgumentParser(description=desc, epilog=epi,
                                 formatter_class=CustomFormatter)
parser.add_argument('outdir', type=str, help='Output directory')
parser.add_argument('--start-date', type=str, default=start_date.strftime('%Y-%m-%d'),
                    help='Minimum date to search for papers')
parser.add_argument('--end-date', type=str, default=datetime.now().strftime('%Y-%m-%d'),
                    help='Maximum date to search for papers')
parser.add_argument('--search-terms', type=str, default=None, nargs='+',
                    help='Optional search term to filter results')
parser.add_argument('--email', type=str, default=None,
                    help='Email address for NCBI API')

# functions
def construct_query(start_date: str, end_date: str, search_terms: list = None) -> str:
    """
    Create an Entrez query string with date range and optional search term.
    Args:
        start_date: Start date (YYYY-MM-DD)
        end_date: End date (YYYY-MM-DD)
        search_term: Optional search term
    Returns:
        Entrez query string
    """
    start = datetime.strptime(start_date, "%Y-%m-%d")
    end = datetime.strptime(end_date, "%Y-%m-%d")
    date_range = f"{start.strftime('%Y/%m/%d')}:{end.strftime('%Y/%m/%d')}[PDAT]"

    if search_terms is None:
        search_terms = [
            'single cell RNA sequencing', 
            'single cell RNA-seq', 
            'single cell transcriptomics'
        ]
    
    if search_terms:
        # add quotes to search terms
        search_terms = [f'"{term}"' for term in search_terms]
        search_terms = ' OR '.join(search_terms)
        query = f'({search_terms}) AND {date_range}'
    else:
        query = date_range
    query += ' AND "human"[Organism]'
    return query

def get_ids(db: DatabaseType, query: str, retmax: int) -> List[str]:
    """
    Search database and return matching IDs.
    Args:
        db: Database to search ("sra" or "gds")
        query: Entrez query string
        retmax: Maximum number of results to return
    Returns:
        List of matching IDs
    """
    try:
        search_handle = Entrez.esearch(db=db, term=query, retmax=retmax)
        search_results = Entrez.read(search_handle)
        search_handle.close()
        return search_results["IdList"]
    except Exception as e:
        print(f"Error searching {db}: {str(e)}")
        return []

def xml_to_dict(element) -> Dict[str, Any]:
    """
    Convert an XML element to a dictionary.
    Args:
        element: XML element
    Returns:
        Dictionary representation of the XML element
    """
    node = {}
    # Include attributes
    if element.attrib:
        node.update(element.attrib)
    # Include child elements
    children = list(element)
    if children:
        child_dict = {}
        for child in children:
            child_result = xml_to_dict(child)
            if child.tag in child_dict:
                # If tag already exists, convert to a list
                if not isinstance(child_dict[child.tag], list):
                    child_dict[child.tag] = [child_dict[child.tag]]
                child_dict[child.tag].append(child_result)
            else:
                child_dict[child.tag] = child_result
        node.update(child_dict)
    else:
        # Include text content
        if element.text and element.text.strip():
            node['text'] = element.text.strip()
    return node

def fetch_dataset(db: DatabaseType, dataset_id: str) -> list:
    """
    Fetch and parse dataset information.
    Args:
        db: Database to search ("sra" or "gds")
        dataset_id: Dataset ID
    Returns:
        Dictionary containing parsed dataset
    """
    # Fetch dataset record
    handle = Entrez.efetch(db=db, id=dataset_id, retmode="xml")
    record = handle.read()
    handle.close()

    # convert to dictionary, if not a string
    if isinstance(record, str):
        return [record]
    root = ET.fromstring(record.decode('utf-8'))
    data = list()
    for child in root:
        data.append(xml_to_dict(child))
    return data

def search_datasets(email: str, start_date: str, end_date: str, db: DatabaseType,
                    search_terms: list = None, retmax: int = 20) -> List[Dict]:
    """
    Search for recent datasets in specified database.
    Args:
        email: Your email address (required by NCBI)
        start_date: Start date (YYYY-MM-DD)
        end_date: End date (YYYY-MM-DD)
        db: Database to search ("sra" or "gds")
        search_terms: Optional search term to filter results
        retmax: Maximum number of results to return
    Returns:
        List of dictionaries containing dataset information
    """
    Entrez.email = email
    
    query = construct_query(start_date, end_date, search_terms)
    dataset_ids = get_ids(db, query, retmax)

    print(f"Found {len(dataset_ids)} datasets.", file=sys.stderr)
    
    datasets = {}
    for id in dataset_ids:
        print(f"  Fetching dataset \"{id}\"...", file=sys.stderr)
        time.sleep(0.5)
        ret = fetch_dataset(db, id)
        if ret:
            datasets[id] = ret

    return datasets

def write_dataset_info(id: str, dataset: Dict, outdir: str) -> None:
    """
    Write dataset information as JSON to a file.
    Args:
        dataset: Dictionary containing dataset information
        outdir: Output directory
    """
    filename = os.path.join(outdir, f"{id}.json")
    with open(filename, "w") as f:
        json.dump(dataset, f, indent=4)
    print(f"Dataset info written to {filename}")

def main(args):
    # status
    print(f"Searching for datasets between {args.start_date} and {args.end_date}...", file=sys.stderr)

    # output
    os.makedirs(args.outdir, exist_ok=True)

    # Search SRA
    print("Searching SRA...", file=sys.stderr)
    sra_datasets = search_datasets(
        args.email, 
        start_date=args.start_date, end_date=args.end_date, 
        db="sra", search_terms=args.search_terms
    )
    for id,dataset in sra_datasets.items():
        write_dataset_info(id, dataset, args.outdir)
    
    # Search GEO
    print("Searching GEO...", file=sys.stderr)
    geo_datasets = search_datasets(
        args.email, 
        start_date=args.start_date, end_date=args.end_date, 
        db="gds", search_terms=args.search_terms
    )
    for id,dataset in geo_datasets.items():
        write_dataset_info(id, dataset, args.outdir)

# Example usage
if __name__ == "__main__":
    args = parser.parse_args()
    main(args)