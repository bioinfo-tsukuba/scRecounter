# ProcessTracker統合ガイド

このドキュメントでは、scRecounterワークフローにProcessTrackerを統合する方法について説明します。

## 概要

ProcessTrackerは、scRecounterの実行状況をローカルPostgreSQLデータベースに記録するシステムです。**1つのExperimental IDにつき1レコード**で、各実験の開始・終了時刻、ステータス、エラーメッセージなどを追跡できます。

### 記録単位
- **1つのExperimental ID** = **1つのscRecounter処理** = **1レコード**
- process_typeは常に`"scRecounter"`
- エラー処理では直近に出力されたエラーメッセージを格納

## 現在の実装状況

ProcessTracker統合は**コメントアウトされた状態**で実装されています。実際に使用するには、該当箇所のコメントを外す必要があります。

## 有効化手順

### 1. ローカルPostgreSQLデータベースのセットアップ

```bash
# PostgreSQLをインストール（Ubuntu/Debian）
sudo apt update
sudo apt install postgresql postgresql-contrib

# データベースとユーザーを作成
sudo -u postgres createdb experimentprocess
sudo -u postgres createuser cellio
sudo -u postgres psql -c "ALTER USER cellio PASSWORD 'cEllIo_process';"
sudo -u postgres psql -c "GRANT ALL PRIVILEGES ON DATABASE experimentprocess TO cellio;"
```

### 2. 環境設定

`.env.local`ファイルが既に存在します：

```bash
LOCAL_DB_HOST=localhost
LOCAL_DB_NAME=experimentprocess
LOCAL_DB_USER=cellio
LOCAL_DB_PASSWORD=cEllIo_process
LOCAL_DB_PORT=5432
```

### 3. ProcessTracker統合の有効化

以下のファイルでコメントアウトされたコードを有効化してください：

#### main.nf
```groovy
// 以下のコメントを外す
experiment_id = "${workflow.runName}_${workflow.start.format('yyyyMMdd_HHmmss')}"
process_type = "scRecounter"
process_id = "version_0.1"

// Start process tracking
"""
python3 ${projectDir}/bin/process_tracker_start.py \\
    --experiment_id ${experiment_id} \\
    --process_type ${process_type} \\
    --process_id ${process_id}
"""

// workflow.onCompleteセクションのコメントも外す
```

**注意**: ワークフローファイル（workflows/）にはサブプロセストラッキングのコードは含まれていません。メインワークフローのみでトラッキングを行います。

## 使用方法

### ProcessTrackerの直接使用

```python
from process_tracker import ProcessTracker

# 初期化（自動的にローカルDBに接続）
tracker = ProcessTracker()

# プロセス開始
process_id = tracker.start_process(
    experiment_id="test_exp_001",
    process_type="scRecounter",
    process_id="version_0.1"
)

# プロセス完了
tracker.finish_process(
    experiment_id="test_exp_001",
    process_type="scRecounter", 
    process_id="version_0.1",
    status=0  # 0=成功, 1=エラー
)
```

### コマンドライン使用

```bash
# プロセス開始
python3 bin/process_tracker_start.py \\
    --experiment_id "test_exp_001" \\
    --process_type "scRecounter" \\
    --process_id "version_0.1"

# プロセス完了
python3 bin/process_tracker_finish.py \\
    --experiment_id "test_exp_001" \\
    --process_type "scRecounter" \\
    --process_id "version_0.1" \\
    --status 0
```

## データベーススキーマ

ProcessTrackerは以下のテーブル構造を使用します：

```sql
CREATE TABLE IF NOT EXISTS experiment_process (
    id INTEGER PRIMARY KEY,
    experiment_id VARCHAR(50) NOT NULL,
    srx_accession VARCHAR(50),
    organism VARCHAR(50),
    analysis_date DATE,
    process_type VARCHAR(50) NOT NULL,
    status INTEGER DEFAULT NULL,  -- NULL=実行中, 0=成功, 1=エラー
    start_datetime TIMESTAMP,
    finish_datetime TIMESTAMP,
    path VARCHAR(500),
    process_id VARCHAR(100),
    error_message TEXT
);
```

## 監視とクエリ

### 実行中プロセスの確認

```python
tracker = ProcessTracker()

# 実行中プロセス一覧
running = tracker.get_pending_processes()
print(running)

# 失敗したプロセス一覧  
failed = tracker.get_failed_processes()
print(failed)

# 成功したプロセス一覧
successful = tracker.get_successful_processes()
print(successful)
```

## トラブルシューティング

### データベース接続エラー
- PostgreSQLサービスが起動していることを確認
- `.env.local`の設定値が正しいことを確認
- データベースとユーザーが作成されていることを確認

### テーブル作成エラー
- ProcessTrackerは初回実行時に自動的にテーブルを作成します
- 権限エラーが発生した場合は、ユーザー権限を確認してください

## エラー処理の特徴

- **包括的なエラー収集**: `workflow.errorMessage`と`workflow.errorReport`の両方を収集
- **直近エラーの保存**: 各Experimental IDに対して最も直近のエラーメッセージを格納
- **エラー情報の統合**: 複数のエラーソースを`|`で結合して保存

## 注意事項

- 現在の実装はローカルPostgreSQLのみをサポートしています
- GCP Cloud SQL関連のコードはコメントアウトされています  
- ProcessTracker統合はオプションです。コメントアウトしたままでもワークフローは正常に動作します
- **1つのExperimental IDにつき1レコード**のみ作成されます（サブプロセストラッキングなし）