#!/bin/bash

# Exit immediately if any command fails
set -e

# Check if GOOGLE_APPLICATION_CREDENTIALS is set
if [[ -z "${GOOGLE_APPLICATION_CREDENTIALS}" ]]; then
  echo "Error: GOOGLE_APPLICATION_CREDENTIALS environment variable is not set."
  exit 1
fi

# Check if the credentials file exists
if [[ ! -f "${GOOGLE_APPLICATION_CREDENTIALS}" ]]; then
  echo "Error: Credentials file '${GOOGLE_APPLICATION_CREDENTIALS}' does not exist."
  exit 1
fi

# Start the Cloud SQL Proxy
echo "Starting Cloud SQL Proxy..."
cloud_sql_proxy --unix-socket=/cloudsql --credentials-file="${GOOGLE_APPLICATION_CREDENTIALS}" c-tc-429521:us-east1:sragent &

# Wait for the Cloud SQL Proxy to initialize
sleep 2

# Run the Python script
/usr/bin/auto-run.py "$@"