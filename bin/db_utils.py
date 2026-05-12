# import
## batteries
import os
import logging
import warnings
from typing import List
## 3rd party
import numpy as np
import psycopg2
import pandas as pd
from psycopg2.extras import execute_values
from psycopg2.extensions import connection
from tempfile import NamedTemporaryFile

# Suppress notifications
warnings.filterwarnings("ignore", message="pandas only supports SQLAlchemy connectable")
logging.getLogger("psycopg2").setLevel(logging.CRITICAL)
# GCP-related logging suppression - COMMENTED OUT for local development
# logging.getLogger("google.auth.transport.requests").setLevel(logging.CRITICAL)
# logging.getLogger("urllib3").setLevel(logging.CRITICAL)
# logging.getLogger("google.auth").setLevel(logging.CRITICAL)

# functions
def db_connect_local() -> connection:
    """Connect to the local PostgreSQL database using environment variables.

    Reads connection parameters from the environment, falling back to defaults
    when variables are absent.

    Returns
    -------
    psycopg2.extensions.connection
        An open database connection.

    Raises
    ------
    Exception
        If the connection cannot be established.
    """
    db_params = {
        'host': os.environ.get("LOCAL_DB_HOST", "localhost"),
        'database': os.environ.get("LOCAL_DB_NAME", "screcounter"),
        'user': os.environ.get("LOCAL_DB_USER", "postgres"),
        'password': os.environ.get("LOCAL_DB_PASSWORD", ""),
        'port': os.environ.get("LOCAL_DB_PORT", "5432"),
        'sslmode': 'disable'
    }
    try:
        conn = psycopg2.connect(**db_params)
        return conn
    except psycopg2.OperationalError as e:
        logging.error(f"Failed to connect to local database: {e}")
        raise Exception(f"Database connection failed: {e}")

# GCP Cloud SQL connection - COMMENTED OUT for local development
# def db_connect() -> connection:
#     """
#     Connect to the sql database using SSL certificates.
#     """
#     # get certs
#     certs = get_db_certs()
#     # connect
#     db_params = {
#         'host': os.environ["GCP_SQL_DB_HOST"],
#         'database': os.environ["GCP_SQL_DB_NAME"],
#         'user': os.environ["GCP_SQL_DB_USERNAME"],
#         'password': os.getenv("GCP_SQL_DB_PASSWORD", get_secret("GCP_SQL_DB_PASSWORD")),
#         'sslmode': 'verify-ca',
#         'sslrootcert': certs["server-ca.pem"],
#         'sslcert': certs["client-cert.pem"],
#         'sslkey': certs["client-key.pem"],
#         'port': '5432',
#         'connect_timeout': 30
#     }
#     conn = psycopg2.connect(**db_params)
#     # delete certs
#     for cert in certs.values():
#         os.remove(cert)
#     return conn

def add_to_log(
        df, sample: str, accession: str, process: str, step: str, status: str, msg: str
        ) -> pd.DataFrame:
    """Append a log entry to an in-memory log DataFrame.

    Truncates messages longer than 200 characters.

    Parameters
    ----------
    df : pd.DataFrame
        Existing log DataFrame to append to.
    sample : str
        Sample name.
    accession : str
        SRA accession identifier.
    process : str
        Name of the process being logged.
    step : str
        Name of the step within the process.
    status : str
        Status string (e.g. 'Success', 'Failure').
    msg : str
        Log message; truncated to 200 characters if longer.

    Returns
    -------
    pd.DataFrame
        The updated log DataFrame.
    """
    if len(msg) > 200:
        msg = str(msg[:(200-3)]) + '...'
    df.loc[len(df)] = [sample, accession, process, step, status, msg]

def sanitize_int_columns(df, min_int=-2**30, max_int=2**30 - 1) -> pd.DataFrame:
    """Cast integer columns to float and replace out-of-range values with NaN.

    PostgreSQL integer bounds differ from Python's; this prevents overflow
    errors when inserting large integer values.

    Parameters
    ----------
    df : pd.DataFrame
        Input DataFrame.
    min_int : int, optional
        Minimum acceptable integer value, by default -(2**30).
    max_int : int, optional
        Maximum acceptable integer value, by default 2**30 - 1.

    Returns
    -------
    pd.DataFrame
        DataFrame with integer columns cast to float and out-of-range values
        replaced by NaN.
    """
    int_cols = df.select_dtypes(include=["int", "int32", "int64"]).columns

    # Cast these columns to float so they can hold NaN values
    for col in int_cols:
        df[col] = df[col].astype(float)

        # Replace out-of-range values with NaN
        df.loc[df[col] < min_int, col] = np.nan
        df.loc[df[col] > max_int, col] = np.nan

    return df

def db_upsert(df: pd.DataFrame, table_name: str, conn: connection) -> None:
    """Insert rows from a DataFrame into a PostgreSQL table, ignoring conflicts.

    Only columns present in both the DataFrame and the target table are
    inserted. Duplicate rows (based on the table's unique constraints) are
    silently skipped via ON CONFLICT DO NOTHING.

    Parameters
    ----------
    df : pd.DataFrame
        Data to insert.
    table_name : str
        Name of the target PostgreSQL table.
    conn : psycopg2.extensions.connection
        Active database connection.

    Raises
    ------
    Exception
        If the INSERT statement fails for any reason other than a conflict.
    """
    # if df is empty, return
    if df.empty:
        return
    # if df is not dataframe, try to convert
    if not isinstance(df, pd.DataFrame):
        try:
            df = pd.DataFrame(df)
        except Exception as e:
            raise Exception(f"Error converting input to DataFrame: {str(e)}")

    # filter to overlapping target columns
    table_columns = get_table_columns(table_name, conn)
    overlapping_columns = list(set(table_columns).intersection(df.columns))
    df = df[overlapping_columns]

    # Sanitize integer columns
    df = sanitize_int_columns(df.copy())

    # Get DataFrame columns
    columns = list(df.columns)

    # Remove duplicates within the DataFrame
    df = df.drop_duplicates(keep='first').copy()

    # Convert DataFrame to list of tuples
    values = [tuple(x) for x in df.to_numpy()]

    # Create simple INSERT statement with DO NOTHING on conflict
    insert_stmt = f"INSERT INTO {table_name} ({', '.join(columns)})"
    insert_stmt += f"\nVALUES %s"
    insert_stmt += f"\nON CONFLICT DO NOTHING"

    # Execute the query
    try:
        with conn.cursor() as cur:
            execute_values(cur, insert_stmt, values)
            conn.commit()
    except Exception as e:
        conn.rollback()
        raise Exception(f"Error uploading data to {table_name}: {str(e)}")

def db_update(df: pd.DataFrame, table_name: str, conn: connection) -> None:
    """Update existing rows in a PostgreSQL table using unique constraint columns as keys.

    Columns participating in the table's unique constraint are used as the
    WHERE condition; all other columns in the DataFrame are updated.

    Parameters
    ----------
    df : pd.DataFrame
        Data containing updated values. Must include the unique-constraint
        columns so that the correct rows can be identified.
    table_name : str
        Name of the target PostgreSQL table.
    conn : psycopg2.extensions.connection
        Active database connection.

    Raises
    ------
    Exception
        If the UPDATE statement fails.
    """
    if df.empty:
        return
    if not isinstance(df, pd.DataFrame):
        df = pd.DataFrame(df)

    # Filter columns to only those that exist in the table
    df = df[list(set(get_table_columns(table_name, conn)).intersection(df.columns))]

    # Sanitize integers, drop duplicates, etc.
    df = sanitize_int_columns(df.copy())
    unique_columns = get_unique_columns(table_name, conn)

    # Get non-unique columns
    columns = list(df.columns)
    non_unique_cols = [c for c in columns if c not in unique_columns]

    # Nothing to update if there are no non-unique columns or no rows
    if not non_unique_cols or df.empty:
        return

    # Convert DataFrame rows to tuples
    values = [tuple(x) for x in df.to_numpy()]

    # Build the WITH data(...) clause
    with_data_cols = ", ".join(columns)
    with_clause = f"WITH data({with_data_cols}) AS (VALUES %s)"

    # Build the UPDATE ... SET ... FROM data ... WHERE ...
    set_clause = ", ".join(f"{col} = data.{col}" for col in non_unique_cols)
    join_condition = " AND ".join(f"t.{uc} = data.{uc}" for uc in unique_columns)

    update_stmt = f"""
    {with_clause}
    UPDATE {table_name} AS t
    SET {set_clause}
    FROM data
    WHERE {join_condition}
    """

    try:
        with conn.cursor() as cur:
            execute_values(cur, update_stmt, values)
        conn.commit()
    except Exception as e:
        conn.rollback()
        raise Exception(f"Error updating data in {table_name}: {str(e)}")

def get_table_columns(table: str, conn: connection) -> List[str]:
    """Return the column names for a table as defined in the database schema.

    Parameters
    ----------
    table : str
        Name of the table to inspect.
    conn : psycopg2.extensions.connection
        Active database connection.

    Returns
    -------
    list of str
        Column names in the order returned by information_schema.columns.
    """
    query = """
    SELECT column_name
    FROM information_schema.columns
    WHERE table_name = %s;
    """

    with conn.cursor() as cur:
        cur.execute(query, (table,))
        columns = cur.fetchall()
    return [col[0] for col in columns]

def get_unique_columns(table: str, conn: connection) -> List[str]:
    """Return the columns forming the most appropriate unique constraint for a table.

    Composite unique constraints are preferred over single-column primary keys.
    Falls back to the primary key if no other constraint is found.

    Parameters
    ----------
    table : str
        Name of the table to inspect.
    conn : psycopg2.extensions.connection
        Active database connection.

    Returns
    -------
    list of str
        Column names that form the selected unique constraint.

    Raises
    ------
    ValueError
        If no unique or primary-key constraint exists on the table.
    """
    query = """
    SELECT c.contype, ARRAY_AGG(a.attname ORDER BY array_position(c.conkey, a.attnum)) as columns
    FROM pg_constraint c
    JOIN pg_class t ON c.conrelid = t.oid
    JOIN pg_attribute a ON a.attrelid = t.oid AND a.attnum = ANY(c.conkey)
    WHERE t.relname = %s
    AND c.contype IN ('p', 'u')  -- primary key or unique constraint
    GROUP BY c.conname, c.contype
    ORDER BY c.contype DESC;  -- 'u'nique before 'p'rimary key
    """

    with conn.cursor() as cur:
        cur.execute(query, (table,))
        constraints = cur.fetchall()

    if not constraints:
        raise ValueError(f"No unique constraints found in table {table}")

    # Prefer composite unique constraints over single-column primary keys
    for constraint_type, columns in constraints:
        if len(columns) > 1 or constraint_type == 'u':
            return columns

    # Fall back to primary key if no other suitable constraint found
    return constraints[0][1]

def get_srx_metadata_limit5(conn):
    """Fetch the first five rows of the srx_metadata table.

    Parameters
    ----------
    conn : psycopg2.extensions.connection
        Active database connection.

    Returns
    -------
    pd.DataFrame
        Up to five rows from srx_metadata.
    """
    query = """
    SELECT * FROM srx_metadata LIMIT 5;
    """
    return pd.read_sql(query, conn)

# main
if __name__ == "__main__":
    from dotenv import load_dotenv
    load_dotenv('.env.local')
    # Test local connection instead of GCP
    # with db_connect_local() as conn:
    #     print(get_srx_metadata_limit5(conn))
