#!/usr/bin/env python
# import
from __future__ import print_function
import os
import re
import sys
import argparse
import logging
from typing import List
## 3rd party
import pandas as pd
import psycopg2
import pandas as pd
from pypika import Query, Table, Criterion
from psycopg2.extras import execute_values
from psycopg2.extensions import connection
## pipeline
from db_utils import db_connect, db_upsert, add_to_log


# logging
logging.basicConfig(format='%(asctime)s - %(message)s', level=logging.DEBUG)

# argparse
class CustomFormatter(argparse.ArgumentDefaultsHelpFormatter,
                      argparse.RawDescriptionHelpFormatter):
    pass

desc = 'Get SRA accessions from the scRecounter database'
epi = """DESCRIPTION:
Get SRA accessions from the scRecounter database. 
Write out the accessions csv table:
- sample: SRX accession
- accession: SRR accession
- organism: organism name
"""
parser = argparse.ArgumentParser(description=desc, epilog=epi,
                                 formatter_class=CustomFormatter)
parser.add_argument('--max-srx', type=int, default=1,
                    help='Max number of srx records to return')
parser.add_argument('--database', type=str, default=["sra", "gds"], nargs="+",
                    help='Only return records from these databases')
parser.add_argument('--outfile', type=str, default="accessions.csv",
                    help='Output file name')

# functions
def db_get_unprocessed_records(
    conn: connection, 
    process: str,
    database: List[str], 
    max_srx: int=3
    ) -> pd.DataFrame:
    """
    Get all suitable unprocessed SRX records, limiting by unique srx_accession values.
    Args:
        conn: Connection to the database.
        database: Name of the database to query.
        max_srx: Maximum number of SRX records to return.
    Returns:
        dataframe of unprocessed SRX records.
    """
    # init tables
    srx_metadata = Table("srx_metadata")
    srx_srr = Table("srx_srr")
    scr_log = Table("screcounter_log")

    # subquery to get srx_accessions
    ## sc-recounter log
    nontarget_srx = Query \
        .from_(scr_log) \
        .select(scr_log.sample) \
        .where(Criterion.all([
            scr_log.process == process,
            scr_log.step == "Final",
            scr_log.status == "Success",
        ])) \
        .distinct()

    ## filtering criteria
    target_srx = Query \
        .from_(srx_metadata) \
        .left_join(nontarget_srx) \
        .on(srx_metadata.srx_accession == nontarget_srx.sample) \
        .select(srx_metadata.srx_accession) \
        .where(Criterion.all([
            nontarget_srx.sample.isnull(),
            srx_metadata.database.isin(database),
            srx_metadata.is_illumina == "yes",
            srx_metadata.is_single_cell == "yes",
            srx_metadata.is_paired_end == "yes",
            ~srx_metadata.tech_10x.isin(["other", "not_applicable"])
        ])) \
        .distinct() \
        .limit(max_srx)

    # main query to obtain the SRR for each SRX and then format the output
    stmt = (
        Query
        .from_(srx_metadata)
        .inner_join(srx_srr)
        .on(srx_metadata.srx_accession == srx_srr.srx_accession)
        .where(
            srx_metadata.srx_accession.isin(target_srx)
        )
        .select(
            srx_metadata.srx_accession.as_("sample"),
            srx_srr.srr_accession.as_("accession"),
            srx_metadata.organism.as_("organism"),
            srx_metadata.tech_10x.as_("tech_10x"),
        )
    )
        
    # fetch as pandas dataframe
    return pd.read_sql(str(stmt), conn)

def main(args):
    # set process name
    process = "Get db accessions"

    # get unprocessed records
    with db_connect() as conn:
        df = db_get_unprocessed_records(conn, process, args.database, max_srx=args.max_srx)
    ## write out records
    df.to_csv(args.outfile, index=False)
    
    # write to log table in scRecounter database
    ## convert df
    df["process"] = process
    df["step"] = "Final"
    df["status"] = "Success"
    df["message"] = "Obtained database accession for processing"

    ## filter columns
    df = df[["sample", "accession", "process", "step", "status", "message"]]
    
    ## upsert log to database
    with db_connect() as conn:
        db_upsert(df, "screcounter_log", conn)

## script main
if __name__ == '__main__':
    args = parser.parse_args()
    main(args)

