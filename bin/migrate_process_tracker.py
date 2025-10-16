#!/usr/bin/env python3
"""
Database schema migration script
Add new columns to experiment_process table
"""

import logging
import sys
from db_utils import db_connect

def migrate_database():
    """Add new columns to existing experiment_process table"""
    
    conn = db_connect()
    
    try:
        with conn.cursor() as cur:
            # Add new columns (only if they don't exist)
            migration_sql = """
            -- Add sample_id column
            DO $$ 
            BEGIN 
                IF NOT EXISTS (SELECT 1 FROM information_schema.columns 
                              WHERE table_name='experiment_process' AND column_name='sample_id') THEN
                    ALTER TABLE experiment_process ADD COLUMN sample_id VARCHAR(200);
                END IF;
            END $$;
            
            -- Add srx_accession column
            DO $$ 
            BEGIN 
                IF NOT EXISTS (SELECT 1 FROM information_schema.columns 
                              WHERE table_name='experiment_process' AND column_name='srx_accession') THEN
                    ALTER TABLE experiment_process ADD COLUMN srx_accession VARCHAR(50);
                END IF;
            END $$;
            
            -- Add organism column
            DO $$ 
            BEGIN 
                IF NOT EXISTS (SELECT 1 FROM information_schema.columns 
                              WHERE table_name='experiment_process' AND column_name='organism') THEN
                    ALTER TABLE experiment_process ADD COLUMN organism VARCHAR(50);
                END IF;
            END $$;
            
            -- Add analysis_date column
            DO $$ 
            BEGIN 
                IF NOT EXISTS (SELECT 1 FROM information_schema.columns 
                              WHERE table_name='experiment_process' AND column_name='analysis_date') THEN
                    ALTER TABLE experiment_process ADD COLUMN analysis_date DATE;
                END IF;
            END $$;
            
            -- Add new indexes
            CREATE INDEX IF NOT EXISTS idx_experiment_process_sample_id ON experiment_process(sample_id);
            CREATE INDEX IF NOT EXISTS idx_experiment_process_srx_accession ON experiment_process(srx_accession);
            """
            
            cur.execute(migration_sql)
            conn.commit()
            
            logging.info("Database migration completed successfully")
            print("✅ Database schema update completed successfully")
            
    except Exception as e:
        conn.rollback()
        logging.error(f"Migration failed: {e}")
        print(f"❌ Migration failed: {e}")
        sys.exit(1)
        
    finally:
        conn.close()

if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    migrate_database()