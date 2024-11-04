#!/usr/bin/env python
# import
import os
import io
import csv
import sys
import argparse
from time import sleep
from typing import List, Dict
from urllib.error import HTTPError
from dotenv import load_dotenv
import pandas as pd
from Bio import Entrez
from pysradb.sraweb import SRAweb


# argparse
class CustomFormatter(argparse.ArgumentDefaultsHelpFormatter,
                      argparse.RawDescriptionHelpFormatter):
    pass

desc = 'Convert accessions to SRR accessions'
epi = """DESCRIPTION:
Convert SRP, GSE, or other accessions to SRR accessions.
If NCBI_API_KEY is set in the environment, it will be used as the API key.
"""
parser = argparse.ArgumentParser(description=desc, epilog=epi,
                                 formatter_class=CustomFormatter)
parser.add_argument('accession_file', type=str, 
                    help='Text file with accessions; 1 per line')
parser.add_argument('--email', type=str, default=None,
                    help='Email address for Entrez')
parser.add_argument('--batch-size', type=int, default=50,
                    help='Batch size for fetching')
parser.add_argument('--outfile', type=str, default='srr_accessions.csv',
                    help='Output file name')

# functions
def load_accessions(accession_file: str) -> List[str]:
    """
    Load accessions from file
    Args:
        accession_file: File with accessions
    Returns:
        List of accessions
    """
    accessions = []
    with open(accession_file) as inF:
        for line in inF:
            line = line.strip().split(',')[0]
            if line == "" or line.startswith("#"):
                continue
            accessions.append(line)
    return accessions

def esearch_batch(db, accession, batch_size = 50, ntries=3, sleep_time=5) -> List[str]:
    """
    Entrez esearch in batches
    Args:
        db: Database to search
        accession: Accession to search
        batch_size: Batch size for fetching
        ntries: Number of tries before giving up
        sleep_time: Sleep time between retries
    Returns:
        List of unique IDs
    """
    print(f"esearch of {db} for: {accession}", file=sys.stderr)
    results = []
    # Initial search to get the total count of records
    handle = Entrez.esearch(db=db, term=accession, usehistory="y", retmax=1)
    record = Entrez.read(handle)
    handle.close()
    results += record["IdList"] if record["IdList"] else None
    total_records = int(record["Count"])
    print(f"  Total records: {total_records}", file=sys.stderr)

    # Retrieve results in batches
    for start in range(0, total_records, batch_size):
        print(f"  Fetching records {start+1}-{min(start+batch_size, total_records)}", file=sys.stderr)
        for i in range(ntries):
            try:
                handle = Entrez.esearch(db=db, term=accession, retstart=start, retmax=batch_size, usehistory="y")
                record = Entrez.read(handle)
                handle.close()
                if "IdList" in record:
                    results += record["IdList"] 
                sleep(0.5)  # comply with NCBI rate limits
            except Exception as e:
                print(f"  Attempt {i+1}/{ntries}: Error encountered: {e}", file=sys.stderr)
                sleep(5 * (i+1))

    # Return unique IDs
    return list(set(results))

def efetch_batch(db, idlist, batch_size=20, rettype="runinfo", retmode="text", ntries=3, sleep_time=5
                 ) -> List[pd.DataFrame]:
    """
    Entrez efetch in batches
    Args:
        db: Database to search
        idlist: List of IDs to fetch
        batch_size: Batch size for fetching
        rettype: Return type
        retmode: Return mode
        ntries: Number of tries before giving up
        sleep_time: Sleep time between retries
    Returns:
        List of dataframes
    """
    print(f"efetch of {db} for: {len(idlist)} IDs", file=sys.stderr)
    results = []
    for start in range(0, len(idlist), batch_size):
        print(f"  Fetching batch {start+1}-{min(start+batch_size, len(idlist))}", file=sys.stderr)
        batch_ids = ",".join(idlist[start:start + batch_size])  # Get current batch of IDs
        batch_result = None
        for i in range(ntries):  # Retry logic for each batch
            try:
                handle = Entrez.efetch(db=db, id=batch_ids, rettype=rettype, retmode=retmode)
                batch_result = handle.read()
                handle.close()
                # convert to dataframe
                df = pd.read_csv(io.StringIO(batch_result.decode('utf-8')))
                results.append(df)
                sleep(0.5)  # comply with NCBI rate limits
                break  # Exit retry loop on success
            except HTTPError as e:
                print(f"  Attempt {i+1}/{ntries}: HTTPError for batch {start}-{start+batch_size}: {e}", file=sys.stderr)
                sleep(sleep_time * (i + 1))  # Progressive wait time before retry
                continue
        if batch_result is None:
            print(f"  Failed to fetch batch {start}-{start+batch_size}", file=sys.stderr)
    return results

def fetch_srr_from_srp(accession, batch_size=50, ntries=3, sleep_time=5) -> pd.DataFrame:
    """
    Fetch SRR accessions from SRP
    Args:
        accession: SRP accession
        batch_size: Batch size for fetching
        ntries: Number of tries before giving up
        sleep_time: Sleep time between retries
    Returns:
        Dataframe with SRR accessions
    """
    # Search the SRA database for the SRP accession
    idlist = esearch_batch("sra", accession, batch_size=batch_size, ntries=ntries, sleep_time=sleep_time)
    # get IDs from record
    if len(idlist) == 0:
        print(f"No records found for accession: {accession}", file=sys.stderr)
        return []
    # Fetch run info to get SRR accessions
    results = efetch_batch("sra", idlist, batch_size=batch_size, ntries=ntries, sleep_time=sleep_time)
    # concat dataframes
    df = pd.concat(results)
    # return specific columns
    to_keep = [
        "Sample", "Run", "Experiment",  "SRAStudy", "BioProject", 
        "spots", "spots_with_mates", "avgLength", "size_MB"
    ]
    df = df[to_keep].rename(columns={
        "Sample" : "sample",
        "Run" : "accession",
        "Experiment" : "experiment",
        "SRAStudy" : "sra_study",
        "BioProject" : "bioproject",
        "avgLength" : "avg_length",
        "size_MB" : "size_mb"
    })
    # getting just unique for for "accession"
    return df.drop_duplicates(subset=["accession"])

def gse_to_srp(accession: str) -> str:
    """
    Use pysradb to convert GSE to SRP
    Args:
        accession: GSE accession
    Returns:   
        SRP accession
    """
    sradb = SRAweb()
    df = sradb.gse_to_srp(
        [accession],
        detailed=False,
        sample_attribute=False,
        expand_sample_attributes=False,
    )
    srp_accession = df["study_accession"].tolist()[0]
    print(f"Converted GSE to SRP: {srp_accession}", file=sys.stderr)
    return srp_accession

def gsm_to_srp(accession: str) -> str:
    """
    Use pysradb to convert GSM to SRP
    Args:
        accession: GSM accession
    Returns:
        SRP accession
    """
    sradb = SRAweb()
    df = sradb.gsm_to_srp(
        [accession],
        detailed=False,
        sample_attribute=False,
        expand_sample_attributes=False,
    )
    srp_accession = df["study_accession"].tolist()[0]
    print(f"Converted GSM to SRP: {srp_accession}", file=sys.stderr)
    return srp_accession

def convert_to_srp(accession: str) -> str:
    """
    Convert GSE or GSM to SRP
    Args:
        accession: GSE or GSM accession
    Returns:
        SRP accession
    """
    if accession.startswith('GSE'):
        try:
            return gse_to_srp(accession)
        except Exception as e:
            print(f"Error converting GSE to SRP: {e}", file=sys.stderr)
            return None
    elif accession.startswith('GSM'):
        try:
            return gsm_to_srp(accession)
        except Exception as e:
            print(f"Error converting GSM to SRP: {e}", file=sys.stderr)
            return None
    else:
        print(f"Accession type not recognized: {accession}", file=sys.stderr)
        return None

def fetch_srr_from_accession(accession: str, batch_size: int) -> List[pd.DataFrame]:
    """
    Fetch SRR accessions from SRP or GSE
    Args:
        accession: SRP or GSE accession
        batch_size: Batch size for fetching
    Returns:
        List of dataframes with SRR accession info
    """
    print(f"#-- Fetching SRR accessions for: {accession} --#", file=sys.stderr)
    if accession.startswith('GSE') or accession.startswith('GSM'):
        # convert GSE to SRP
        srp_accession = convert_to_srp(accession)
        df = fetch_srr_from_srp(srp_accession)
    elif accession.startswith('SRP'):
        # fetch SRR from SRP
        df = fetch_srr_from_srp(accession)
    else:
        print(f"Accession type not recognized: {accession}", file=sys.stderr)
        return None
    # add query accession
    df["query_accession"] = accession
    # move query accession to first column
    cols = df.columns.tolist()
    cols = cols[-1:] + cols[:-1]
    return df[cols]

def main(args):
    # load accessions 
    accessions = load_accessions(args.accession_file)

    # set email
    if args.email:
        Entrez.email = args.email
    # set API key
    if 'NCBI_API_KEY' in os.environ:
        Entrez.api_key = os.environ['NCBI_API_KEY']

    # get SRR accessions
    srr_accessions = []
    for accession in accessions:
        srr_accessions.append(
            fetch_srr_from_accession(accession, batch_size=args.batch_size)
        )

    # concat list of dataframes
    srr_accessions = pd.concat(srr_accessions)

    # write table
    srr_accessions.to_csv(args.outfile, sep=',', index=False)
    print(f"Saved SRR accessions to: {args.outfile}", file=sys.stderr)


## script main
if __name__ == '__main__':
    args = parser.parse_args()
    load_dotenv()
    main(args)