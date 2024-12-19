# import
## batteries
import os
import warnings
from typing import List
## 3rd party
import psycopg2
import pandas as pd
from psycopg2.extras import execute_values
from psycopg2.extensions import connection

# Suppress the specific warning
warnings.filterwarnings("ignore", message="pandas only supports SQLAlchemy connectable")

# functions
def db_connect() -> connection:
    """Connect to the sql database"""
    PSWD_VAR_NAME = f"GCP_SQL_DB_PASSWORD_{os.environ['GCP_SQL_DB_TENANT']}"
    db_params = {
        'host': os.environ["GCP_SQL_DB_HOST"],
        'database': os.environ["GCP_SQL_DB_NAME"],
        'user': os.environ["GCP_SQL_DB_USERNAME"],
        'password': os.getenv(PSWD_VAR_NAME, get_secret(PSWD_VAR_NAME)),
        'port': '5432',
        'connect_timeout': 30
    }
    return psycopg2.connect(**db_params)

def db_upsert(df: pd.DataFrame, table_name: str, conn: connection) -> None:
    """
    Upload a pandas DataFrame to PostgreSQL, performing an upsert operation.
    If records exist (based on unique constraints), update them; otherwise insert new records.
    Args:
        df: pandas DataFrame to upload
        table_name: name of the target table
        conn: psycopg2 connection object
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

    # Get DataFrame columns
    columns = list(df.columns)
    
    # Create ON CONFLICT clause based on unique constraints
    unique_columns = get_unique_columns(table_name, conn)

    # Exclude 'id' column from the upsert
    if "id" in columns:
        df = df.drop(columns=["id"])
        columns.remove("id")

    # Drop duplicate records based on unique columns
    df.drop_duplicates(subset=unique_columns, keep='first', inplace=True)

    # Convert DataFrame to list of tuples
    values = [tuple(x) for x in df.to_numpy()]

    # Create the INSERT statement with ON CONFLICT clause
    insert_stmt = f"INSERT INTO {table_name} ({', '.join(columns)})"
    insert_stmt += f"\nVALUES %s"

    # Add DO UPDATE SET clause for non-unique columns
    do_update_set = [col for col in columns if col not in unique_columns]
    if do_update_set:
        do_update_set = ', '.join(f"{col} = EXCLUDED.{col}" for col in do_update_set)
        insert_stmt += f"\nON CONFLICT ({', '.join(unique_columns)})"
        insert_stmt += f"\nDO UPDATE SET {do_update_set}"
    else:
        # if no non-unique columns, add DO NOTHING clause
        insert_stmt += f"\nON CONFLICT ({', '.join(unique_columns)}) DO NOTHING"

    # Execute the query
    try:
        with conn.cursor() as cur:
            execute_values(cur, insert_stmt, values)
            conn.commit()
    except Exception as e:
        conn.rollback()
        raise Exception(f"Error uploading data to {table_name}: {str(e)}")


def get_unique_columns(table: str, conn: connection) -> List[str]:
    """
    Get all unique constraint columns for a table from the database schema.
    Prioritizes composite unique constraints over primary keys.
    
    Args:
        table: Name of the table
        conn: Database connection
        
    Returns:
        List of column names that form the most appropriate unique constraint
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

def get_secret(secret_id: str) -> str:
    """
    Fetch secret from GCP Secret Manager.
    Rquired environment variables: GCP_PROJECT_ID, GOOGLE_APPLICATION_CREDENTIALS
    Args:
        secret_id: The secret id
    Returns:
        The secret value
    """
    from google.auth import default
    from google.cloud import secretmanager

    _, project_id = default()  # Use default credentials; project_id is inferred
    if not project_id:
        project_id = os.environ["GCP_PROJECT_ID"]
    name = f"projects/{project_id}/secrets/{secret_id}/versions/latest"
    client = secretmanager.SecretManagerServiceClient()
    response = client.access_secret_version(request={"name": name})
    return response.payload.data.decode('UTF-8')