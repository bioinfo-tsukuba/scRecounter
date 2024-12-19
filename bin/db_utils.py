# import
## batteries
import os
import warnings
from typing import List, Dict, Any, Tuple, Optional
## 3rd party
import psycopg2
import pandas as pd
from pypika import Query, Table, Field, Column, Criterion
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