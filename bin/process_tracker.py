# import
## batteries
import hashlib
import json
import logging
import time
from datetime import datetime, timedelta
from typing import Optional, Dict, Any, List, Union
## 3rd party
import pandas as pd
from psycopg2.extensions import connection
## local
from db_utils import db_connect, db_upsert, db_update

class ProcessTracker:
    """scRecounter プロセス実行管理クラス"""
    
    def __init__(self, conn: Optional[connection] = None):
        self.conn = conn or db_connect()
        self._ensure_table_exists()
    
    def _ensure_table_exists(self):
        """テーブルが存在しない場合は作成"""
        create_table_sql = """
        CREATE TABLE IF NOT EXISTS experiment_process (
            id INTEGER PRIMARY KEY,
            experiment_id VARCHAR(50) NOT NULL,
            sample_id VARCHAR(200),
            srx_accession VARCHAR(50),
            organism VARCHAR(50),
            analysis_date DATE,
            process_type VARCHAR(50) NOT NULL,
            status INTEGER DEFAULT NULL,
            start_datetime TIMESTAMP,
            finish_datetime TIMESTAMP,
            path VARCHAR(500),
            process_id VARCHAR(100),
            error_message TEXT
        );
        
        CREATE INDEX IF NOT EXISTS idx_experiment_process_status ON experiment_process(status);
        CREATE INDEX IF NOT EXISTS idx_experiment_process_experiment_id ON experiment_process(experiment_id);
        CREATE INDEX IF NOT EXISTS idx_experiment_process_sample_id ON experiment_process(sample_id);
        CREATE INDEX IF NOT EXISTS idx_experiment_process_srx_accession ON experiment_process(srx_accession);
        CREATE INDEX IF NOT EXISTS idx_experiment_process_process_type ON experiment_process(process_type);
        """
        with self.conn.cursor() as cur:
            cur.execute(create_table_sql)
            self.conn.commit()
    
    def _generate_sample_id(self, srx_accession: str, organism: str) -> str:
        """SRXアクセッションからユニークなサンプルIDを生成"""
        analysis_date = datetime.now().strftime('%Y%m%d')
        timestamp = int(time.time())
        return f"{srx_accession}_{organism}_{analysis_date}_{timestamp}"
    
    def _generate_id(self, experiment_id: str, process_type: str) -> int:
        """Experiment IDとProcess TypeからハッシュIDを生成"""
        hash_str = f"{experiment_id}_{process_type}"
        return int(hashlib.md5(hash_str.encode()).hexdigest()[:8], 16) % (2**31 - 1)
    
    def start_process(self, experiment_id: str, process_type: str = "scRecounter", 
                     process_id: Optional[str] = None, path: Optional[str] = None,
                     srx_accession: Optional[str] = None, organism: Optional[str] = None) -> int:
        """プロセス開始"""
        id = self._generate_id(experiment_id, process_type)
        
        # サンプルIDを生成（新しいパラメータが提供された場合のみ）
        sample_id = None
        if srx_accession and organism:
            sample_id = self._generate_sample_id(srx_accession, organism)
        
        process_data = pd.DataFrame([{
            'id': id,
            'experiment_id': experiment_id,
            'sample_id': sample_id,
            'srx_accession': srx_accession,
            'organism': organism,
            'analysis_date': datetime.now().date() if srx_accession else None,
            'process_type': process_type,
            'process_id': process_id,
            'path': path,
            'start_datetime': datetime.now(),
            'status': None  # 実行中
        }])
        
        db_upsert(process_data, 'experiment_process', self.conn)
        logging.info(f"Started process: {experiment_id} - {process_type}")
        return id
    
    def finish_process(self, experiment_id: str, process_type: str = "scRecounter", 
                      status: int = 0, path: Optional[str] = None, 
                      error_message: Optional[str] = None):
        """プロセス完了"""
        id = self._generate_id(experiment_id, process_type)
        finish_time = datetime.now()
        
        update_data = pd.DataFrame([{
            'id': id,
            'experiment_id': experiment_id,
            'process_type': process_type,
            'status': status,
            'finish_datetime': finish_time,
            'path': path,
            'error_message': error_message
        }])
        
        db_upsert(update_data, 'experiment_process', self.conn)
        status_text = "SUCCESS" if status == 0 else "ERROR"
        logging.info(f"Finished process: {experiment_id} - {process_type} - {status_text}")
    
    def get_process_status(self, experiment_id: str, process_type: str = "scRecounter") -> Dict[str, Any]:
        """プロセス状態取得"""
        query = """
        SELECT * FROM experiment_process 
        WHERE experiment_id = %s AND process_type = %s
        ORDER BY created_at DESC LIMIT 1
        """
        result = pd.read_sql(query, self.conn, params=[experiment_id, process_type])
        return result.to_dict('records')[0] if not result.empty else {}
    
    def get_failed_processes(self) -> pd.DataFrame:
        """エラーで終了したプロセス一覧"""
        query = "SELECT * FROM experiment_process WHERE status = 1"
        return pd.read_sql(query, self.conn)
    
    def get_pending_processes(self) -> pd.DataFrame:
        """未完了プロセス一覧"""
        query = "SELECT * FROM experiment_process WHERE status IS NULL"
        return pd.read_sql(query, self.conn)
    
    def get_successful_processes(self) -> pd.DataFrame:
        """正常完了プロセス一覧"""
        query = "SELECT * FROM experiment_process WHERE status = 0"
        return pd.read_sql(query, self.conn)
    
    def get_all_processes(self, experiment_id: Optional[str] = None) -> pd.DataFrame:
        """全プロセス一覧"""
        if experiment_id:
            query = "SELECT * FROM experiment_process WHERE experiment_id = %s ORDER BY created_at DESC"
            return pd.read_sql(query, self.conn, params=[experiment_id])
        else:
            query = "SELECT * FROM experiment_process ORDER BY created_at DESC"
            return pd.read_sql(query, self.conn)
    
    def get_processes_by_srx(self, srx_accession: str) -> pd.DataFrame:
        """SRXアクセッションで全プロセス取得"""
        query = """
        SELECT * FROM experiment_process 
        WHERE srx_accession = %s 
        ORDER BY created_at DESC
        """
        return pd.read_sql(query, self.conn, params=[srx_accession])
    
    def get_processes_by_sample_id(self, sample_id: str) -> pd.DataFrame:
        """サンプルIDで全プロセス取得"""
        query = """
        SELECT * FROM experiment_process 
        WHERE sample_id = %s 
        ORDER BY created_at DESC
        """
        return pd.read_sql(query, self.conn, params=[sample_id])
    
    def start_batch_processes(self, processes: List[Dict[str, Any]]) -> List[int]:
        """複数プロセスの一括開始"""
        ids = []
        for process_info in processes:
            process_id = self.start_process(**process_info)
            ids.append(process_id)
        return ids
    
    def finish_batch_processes(self, processes: List[Dict[str, Any]]) -> None:
        """複数プロセスの一括完了"""
        for process_info in processes:
            self.finish_process(**process_info)
    
    def get_execution_stats(self, process_type: Optional[str] = None, 
                           days_back: int = 30) -> Dict[str, Any]:
        """実行統計を取得"""
        date_filter = datetime.now() - timedelta(days=days_back)
        
        base_query = """
        SELECT 
            process_type,
            COUNT(*) as total_processes,
            COUNT(CASE WHEN status = 0 THEN 1 END) as successful,
            COUNT(CASE WHEN status = 1 THEN 1 END) as failed,
            COUNT(CASE WHEN status IS NULL THEN 1 END) as pending,
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
        """エラー概要取得"""
        date_filter = datetime.now() - timedelta(days=days_back)
        query = """
        SELECT 
            experiment_id,
            process_type,
            error_message,
            start_datetime,
            finish_datetime,
            execution_time_seconds
        FROM experiment_process 
        WHERE status = 1 AND created_at >= %s
        ORDER BY finish_datetime DESC
        """
        return pd.read_sql(query, self.conn, params=[date_filter])
    
    def retry_failed_process(self, experiment_id: str, process_type: str = "scRecounter") -> int:
        """失敗したプロセスの再実行"""
        # 新しいprocess_idで再実行
        new_process_id = f"{process_type}_{datetime.now().strftime('%Y%m%d_%H%M%S')}_retry"
        return self.start_process(experiment_id, process_type, new_process_id)
    
    def get_running_processes(self) -> pd.DataFrame:
        """実行中プロセス一覧（24時間以上実行中のものを含む）"""
        query = """
        SELECT *,
               EXTRACT(EPOCH FROM (NOW() - start_datetime))/3600 as hours_running
        FROM experiment_process 
        WHERE status IS NULL
        ORDER BY start_datetime DESC
        """
        return pd.read_sql(query, self.conn)
    
    def cleanup_old_records(self, days_to_keep: int = 90) -> int:
        """古いレコードのクリーンアップ"""
        cutoff_date = datetime.now() - timedelta(days=days_to_keep)
        query = "DELETE FROM experiment_process WHERE created_at < %s AND status IS NOT NULL"
        
        with self.conn.cursor() as cur:
            cur.execute(query, [cutoff_date])
            deleted_count = cur.rowcount
            self.conn.commit()
        
        logging.info(f"Cleaned up {deleted_count} old process records")
        return deleted_count
    
    def close(self):
        """データベース接続を閉じる"""
        if self.conn:
            self.conn.close()