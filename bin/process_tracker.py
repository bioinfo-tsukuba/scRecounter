# import
import hashlib
import json
import logging
import time
from datetime import datetime, timedelta
from typing import Optional, Dict, Any, List, Union
import pandas as pd
from psycopg2.extensions import connection
from db_utils import db_connect_local, db_upsert, db_update

class ProcessTracker:
    """Track and manage scRecounter process execution state in the database."""

    def __init__(self, conn: Optional[connection] = None):
        self.conn = conn or db_connect_local()
        # self._ensure_table_exists()

    def _generate_sample_id(self, srx_accession: str, organism: str) -> str:
        """Generate a unique sample ID from SRX accession and organism.

        Parameters
        ----------
        srx_accession : str
            SRX accession number.
        organism : str
            Organism name.

        Returns
        -------
        str
            Unique sample ID combining accession, organism, date, and timestamp.
        """
        analysis_date = datetime.now().strftime('%Y%m%d')
        timestamp = int(time.time())
        return f"{srx_accession}_{organism}_{analysis_date}_{timestamp}"

    def _generate_id(self, experiment_id: str, process_type: str, process_id: str) -> str:
        """Generate a unique record ID from experiment ID, process type, and process ID.

        Parameters
        ----------
        experiment_id : str
            Unique experiment identifier.
        process_type : str
            Type of process (e.g. 'scRecounter').
        process_id : str
            Process version or run identifier.

        Returns
        -------
        str
            Concatenated unique ID string.
        """
        process_str = f"{experiment_id}_{process_type}_{process_id}"
        return process_str

    def check_existing_process(self, experiment_id: str, process_type: str = "scRecounter",
                              process_id: Optional[str] = None) -> bool:
        """Check whether a successfully completed process record already exists.

        Records with status=1 (failed) or status=2 (running/interrupted) are treated
        as eligible for re-execution and will return False.

        Parameters
        ----------
        experiment_id : str
            Unique experiment identifier.
        process_type : str, optional
            Process type label, by default 'scRecounter'.
        process_id : str, optional
            Process version or run identifier.

        Returns
        -------
        bool
            True only if a record with status=0 (success) exists for this ID.
        """
        id = self._generate_id(experiment_id, process_type, process_id)

        query = """
        SELECT COUNT(*) as count FROM experiment_process
        WHERE id = %s AND status = 0
        """
        result = pd.read_sql(query, self.conn, params=[id])
        return result['count'].iloc[0] > 0

    def start_process(self, experiment_id: str, process_type: str = "scRecounter",
                     process_id: Optional[str] = None, path: Optional[str] = None,
                     srx_accession: Optional[str] = None, organism: Optional[str] = None,
                     accessions_file: Optional[str] = None) -> str:
        """Create a new process record or reset a previous incomplete record.

        Uses INSERT ... ON CONFLICT DO UPDATE so that records left in status=2
        (running/interrupted) or status=1 (failed) are reset for re-execution.
        Records with status=0 (success) are never overwritten.

        Parameters
        ----------
        experiment_id : str
            Unique experiment identifier.
        process_type : str, optional
            Process type label, by default 'scRecounter'.
        process_id : str, optional
            Process version or run identifier.
        path : str, optional
            Output path for the process results.
        srx_accession : str, optional
            SRX accession number.
        organism : str, optional
            Organism name.
        accessions_file : str, optional
            Path to the accessions input file.

        Returns
        -------
        str
            The generated record ID.
        """
        id = self._generate_id(experiment_id, process_type, process_id)
        now = datetime.now()

        with self.conn.cursor() as cur:
            cur.execute("""
                INSERT INTO experiment_process
                    (id, experiment_id, srx_accession, organism, analysis_date,
                     process_type, process_id, path, start_datetime, status, accessions_file)
                VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s)
                ON CONFLICT (id) DO UPDATE SET
                    start_datetime  = EXCLUDED.start_datetime,
                    status          = EXCLUDED.status,
                    finish_datetime = NULL
                WHERE experiment_process.status != 0
            """, [id, experiment_id, srx_accession, organism,
                  now.date(), process_type, process_id,
                  path, now, 2, accessions_file])
            self.conn.commit()

        logging.info(f"Started process: {experiment_id} - {process_type} - {process_id}")
        return id

    def finish_process(self, experiment_id: str, process_type: str = "scRecounter",
                      process_id: Optional[str] = None, status: int = 0, path: Optional[str] = None,
                      error_message: Optional[str] = None):
        """Record the completion of a process with its final status.

        Parameters
        ----------
        experiment_id : str
            Unique experiment identifier.
        process_type : str, optional
            Process type label, by default 'scRecounter'.
        process_id : str, optional
            Process version or run identifier. Reads from VERSION file if not provided.
        status : int, optional
            Exit status: 0 for success, 1 for failure, by default 0.
        path : str, optional
            Output path for the process results.
        error_message : str, optional
            Error description if the process failed.
        """
        if process_id is None:
            import os
            try:
                script_dir = os.path.dirname(os.path.abspath(__file__))
                project_dir = os.path.dirname(script_dir)
                version_file = os.path.join(project_dir, 'VERSION')
                with open(version_file, 'r') as f:
                    process_id = f.read().strip()
            except FileNotFoundError:
                process_id = "version_0.1"  # fallback
        id = self._generate_id(experiment_id, process_type, process_id)
        finish_time = datetime.now()

        update_data = pd.DataFrame([{
            'id': id,
            'experiment_id': experiment_id,
            'process_type': process_type,
            'status': status,
            'finish_datetime': finish_time,
            'path': path
        }])

        db_update(update_data, 'experiment_process', self.conn)
        status_text = "SUCCESS" if status == 0 else "ERROR"
        logging.info(f"Finished process: {experiment_id} - {process_type} - {process_id} - {status_text}")

    def get_process_status(self, experiment_id: str, process_id: str, process_type: str = "scRecounter") -> Dict[str, Any]:
        """Retrieve the most recent status record for a process.

        Parameters
        ----------
        experiment_id : str
            Unique experiment identifier.
        process_id : str
            Process version or run identifier.
        process_type : str, optional
            Process type label, by default 'scRecounter'.

        Returns
        -------
        dict
            Row as a dictionary, or empty dict if no record found.
        """
        query = """
        SELECT * FROM experiment_process
        WHERE experiment_id = %s AND process_type = %s　AND process_id = %s
        ORDER BY created_at DESC LIMIT 1
        """
        result = pd.read_sql(query, self.conn, params=[experiment_id, process_type, process_id])
        return result.to_dict('records')[0] if not result.empty else {}

    def get_failed_processes(self) -> pd.DataFrame:
        """Return all process records that ended with an error (status=1).

        Returns
        -------
        pd.DataFrame
            Rows from experiment_process where status = 1.
        """
        query = "SELECT * FROM experiment_process WHERE status = 1"
        return pd.read_sql(query, self.conn)

    def get_pending_processes(self) -> pd.DataFrame:
        """Return all process records currently in running or interrupted state (status=2).

        Returns
        -------
        pd.DataFrame
            Rows from experiment_process where status = 2.
        """
        query = "SELECT * FROM experiment_process WHERE status = 2"
        return pd.read_sql(query, self.conn)

    def get_successful_processes(self) -> pd.DataFrame:
        """Return all process records that completed successfully (status=0).

        Returns
        -------
        pd.DataFrame
            Rows from experiment_process where status = 0.
        """
        query = "SELECT * FROM experiment_process WHERE status = 0"
        return pd.read_sql(query, self.conn)

    def get_all_processes(self, experiment_id: Optional[str] = None) -> pd.DataFrame:
        """Return all process records, optionally filtered by experiment ID.

        Parameters
        ----------
        experiment_id : str, optional
            If provided, restricts results to this experiment.

        Returns
        -------
        pd.DataFrame
            Matching rows from experiment_process ordered by creation time descending.
        """
        if experiment_id:
            query = "SELECT * FROM experiment_process WHERE experiment_id = %s ORDER BY created_at DESC"
            return pd.read_sql(query, self.conn, params=[experiment_id])
        else:
            query = "SELECT * FROM experiment_process ORDER BY created_at DESC"
            return pd.read_sql(query, self.conn)

    def get_processes_by_srx(self, srx_accession: str) -> pd.DataFrame:
        """Return all process records associated with a given SRX accession.

        Parameters
        ----------
        srx_accession : str
            SRX accession number to look up.

        Returns
        -------
        pd.DataFrame
            Matching rows ordered by creation time descending.
        """
        query = """
        SELECT * FROM experiment_process
        WHERE srx_accession = %s
        ORDER BY created_at DESC
        """
        return pd.read_sql(query, self.conn, params=[srx_accession])

    def get_processes_by_sample_id(self, id: str) -> pd.DataFrame:
        """Return all process records for a given sample ID.

        Parameters
        ----------
        id : str
            Sample ID to look up.

        Returns
        -------
        pd.DataFrame
            Matching rows ordered by creation time descending.
        """
        query = """
        SELECT * FROM experiment_process
        WHERE id = %s
        ORDER BY created_at DESC
        """
        return pd.read_sql(query, self.conn, params=[id])

    def start_batch_processes(self, processes: List[Dict[str, Any]]) -> List[str]:
        """Start multiple processes in sequence.

        Parameters
        ----------
        processes : list of dict
            Each dict is passed as keyword arguments to start_process().

        Returns
        -------
        list of str
            Generated record IDs for each started process.
        """
        ids = []
        for process_info in processes:
            process_id = self.start_process(**process_info)
            ids.append(process_id)
        return ids

    def finish_batch_processes(self, processes: List[Dict[str, Any]]) -> None:
        """Finish multiple processes in sequence.

        Parameters
        ----------
        processes : list of dict
            Each dict is passed as keyword arguments to finish_process().
        """
        for process_info in processes:
            self.finish_process(**process_info)

    def get_execution_stats(self, process_type: Optional[str] = None,
                           days_back: int = 30) -> Dict[str, Any]:
        """Retrieve aggregated execution statistics for recent processes.

        Parameters
        ----------
        process_type : str, optional
            If provided, restricts statistics to this process type.
        days_back : int, optional
            Number of days of history to include, by default 30.

        Returns
        -------
        list of dict
            Each dict contains counts (successful, failed, running) and timing
            statistics (avg, max, min execution time) per process_type.
        """
        date_filter = datetime.now() - timedelta(days=days_back)

        base_query = """
        SELECT
            process_type,
            COUNT(*) as total_processes,
            COUNT(CASE WHEN status = 0 THEN 1 END) as successful,
            COUNT(CASE WHEN status = 1 THEN 1 END) as failed,
            COUNT(CASE WHEN status = 2 THEN 1 END) as running,
            AVG(execution_time_seconds) as avg_execution_time,
            MAX(execution_time_seconds) as max_execution_time,
            MIN(execution_time_seconds) as min_execution_time
        FROM experiment_process
        WHERE created_at >= %s
        """

        if process_type:
            query = base_query + " AND process_type = %s GROUP BY process_type"
            params = [date_filter, process_type]
        else:
            query = base_query + " GROUP BY process_type"
            params = [date_filter]

        result = pd.read_sql(query, self.conn, params=params)
        return result.to_dict('records')

    def get_error_summary(self, days_back: int = 7) -> pd.DataFrame:
        """Return a summary of failed processes within a recent time window.

        Parameters
        ----------
        days_back : int, optional
            Number of days of history to include, by default 7.

        Returns
        -------
        pd.DataFrame
            Failed process records with error messages and timing, ordered by
            finish time descending.
        """
        date_filter = datetime.now() - timedelta(days=days_back)
        query = """
        SELECT
            experiment_id,
            process_type,
            start_datetime,
            finish_datetime,
            execution_time_seconds
        FROM experiment_process
        WHERE status = 1 AND created_at >= %s
        ORDER BY finish_datetime DESC
        """
        return pd.read_sql(query, self.conn, params=[date_filter])

    def retry_failed_process(self, experiment_id: str, process_type: str = "scRecounter") -> str:
        """Start a new retry attempt for a failed process using a timestamped process ID.

        Parameters
        ----------
        experiment_id : str
            Unique experiment identifier.
        process_type : str, optional
            Process type label, by default 'scRecounter'.

        Returns
        -------
        str
            The generated record ID for the retry attempt.
        """
        new_process_id = f"{process_type}_{datetime.now().strftime('%Y%m%d_%H%M%S')}_retry"
        return self.start_process(experiment_id, process_type, new_process_id)

    def get_running_processes(self) -> pd.DataFrame:
        """Return all process records currently in running state (status=2).

        Includes the number of hours each process has been running.

        Returns
        -------
        pd.DataFrame
            Rows with status=2, augmented with hours_running, ordered by
            start time descending.
        """
        query = """
        SELECT *,
               EXTRACT(EPOCH FROM (NOW() - start_datetime))/3600 as hours_running
        FROM experiment_process
        WHERE status = 2
        ORDER BY start_datetime DESC
        """
        return pd.read_sql(query, self.conn)

    def mark_as_interrupted(self, process_type: str, process_id: str, since: datetime) -> int:
        """Update running records (status=2) to failed (status=1) for a completed workflow run.

        Called from workflow.onComplete to ensure that accessions dropped from the
        Nextflow channel (e.g. due to errorStrategy 'ignore') are recorded as failed
        rather than left indefinitely in the running state.

        Parameters
        ----------
        process_type : str
            Process type label used to scope the update.
        process_id : str
            Process version string used to scope the update.
        since : datetime
            Workflow start time; only records with start_datetime >= this value
            are updated.

        Returns
        -------
        int
            Number of records updated.
        """
        with self.conn.cursor() as cur:
            cur.execute("""
                UPDATE experiment_process
                SET status          = 1,
                    finish_datetime = NOW()
                WHERE status        = 2
                  AND process_type  = %s
                  AND process_id    = %s
                  AND start_datetime >= %s
            """, [process_type, process_id, since])
            updated = cur.rowcount
            self.conn.commit()
        logging.info(f"Marked {updated} processes as interrupted")
        return updated

    def cleanup_old_records(self, days_to_keep: int = 90) -> int:
        """Delete completed process records older than a specified retention period.

        Only records with status=0 (success) or status=1 (failed) are deleted.
        Records with status=2 (running) are retained regardless of age.

        Parameters
        ----------
        days_to_keep : int, optional
            Retention period in days, by default 90.

        Returns
        -------
        int
            Number of records deleted.
        """
        cutoff_date = datetime.now() - timedelta(days=days_to_keep)
        query = "DELETE FROM experiment_process WHERE created_at < %s AND status IN (0, 1)"

        with self.conn.cursor() as cur:
            cur.execute(query, [cutoff_date])
            deleted_count = cur.rowcount
            self.conn.commit()

        logging.info(f"Cleaned up {deleted_count} old process records")
        return deleted_count

    def close(self):
        """Close the database connection."""
        if self.conn:
            self.conn.close()
