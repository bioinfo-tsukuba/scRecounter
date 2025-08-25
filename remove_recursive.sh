#!/usr/bin/bash
#
# 再帰的にディレクトリを探索し、
# ・特定の拡張子 の通常ファイル（シンボリックリンクは除く）を削除
#

set -euo pipefail

# 引数チェック
if [ $# -ne 2 ]; then
    echo "Usage: $0 <TARGET_DIR>"
    exit 1
fi

TARGET_DIR="$1"
format="$2"

# ディレクトリ存在チェック
if [ ! -d "$TARGET_DIR" ]; then
    echo "Error: Directory '$TARGET_DIR' not found." >&2
    exit 1
fi

# find で再帰検索：
# -type f        → 通常ファイルのみ（シンボリックリンクは除外）
# -name '*.fastq' → 拡張子が .fastq のファイル
find "$TARGET_DIR" -type f -name "*.${format}" -print0 \
  | xargs -0 --no-run-if-empty rm

echo "Done: Remove all .fastq files under '$TARGET_DIR' (symlinks ignored)."
