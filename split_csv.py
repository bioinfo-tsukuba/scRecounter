#!/usr/bin/env python3
"""
Minimal splitter for NCBI SRA SraRunInfo.csv

Assumptions (fixed):
- Input CSV has columns: Experiment, Run, size_MB.
- size_MB is in megabytes (MB), 1 GB = 1000 MB.
- Keep all rows for the same Experiment in the same output file.
- Packing strategy: input-order, sequential fill under limit.
- Outputs are written to the current directory as:
    <prefix>_001.csv, <prefix>_002.csv, ... (columns: Experiment, Run)
    summary.csv (per-file totals)

Arguments (exactly 3):
  1) input_csv : path to SraRunInfo.csv
  2) limit_gb  : per-file size limit in GB (float)
  3) prefix    : output filename prefix

Exit behavior:
- If a single Experiment total exceeds the limit, the script fails.
"""

from __future__ import annotations

import argparse
import sys
from dataclasses import dataclass
from typing import List, Dict

import pandas as pd


@dataclass
class Bin:
    experiments: List[str]
    total_mb: float

def pack_inorder(exp_order: List[str],
                 exp_totals: Dict[str, float],
                 limit_mb: float,
                 allow_oversized: bool,) -> List[Bin]:
    bins: List[Bin] = []
    current = Bin(experiments=[], total_mb=0.0)

    # currentをリセット
    def flush_current():
        nonlocal current
        if current.experiments:
            bins.append(current)
            current = Bin(experiments=[], total_mb=0.0)

    # 各Experimentに対して解析
    for exp in exp_order:
        sz = float(exp_totals.get(exp, 0.0))
        if sz > limit_mb:
            if not allow_oversized:
                raise ValueError(
                    f"Experiment '{exp}' total {sz:.2f} MB exceeds limit {limit_mb:.2f} MB"
                )
            # finalize current and place oversized as a single bin
            flush()
            bins.append(Bin(experiments=[exp], total_mb=sz))
            continue

        # 合計ファイルサイズが大きい場合
        if not current.experiments or current.total_mb + sz <= limit_mb:
            current.experiments.append(exp)
            current.total_mb += sz
        # 合計ファイルサイズが小さい場合
        else:
            flush_current()
            current.experiments.append(exp)
            current.total_mb = sz

    # 最後に残ったbinを更新
    flush_current()
    return bins

def write_outputs(df: pd.DataFrame, bins: List[Bin], prefix: str, organism: str) -> pd.DataFrame:
    rows = []
    for i, b in enumerate(bins, start=1):
        subset = (df[df["Experiment"].isin(b.experiments)][["Experiment", "Run"]]
                  .rename(columns={"Experiment": "sample", "Run": "accession"})
                  .assign(organism=str(organism))[["sample", "accession", "organism"]])
        out_path = f"{prefix}_{i:03d}.csv"
        subset.to_csv(out_path, index=False)
        rows.append(
            {
                "file": out_path,
                "experiments": len(b.experiments),
                "runs": int(subset.shape[0]),
                "total_MB": float(b.total_mb),
                "total_GB": float(b.total_mb) / 1000.0,
            }
        )
    summary = pd.DataFrame(rows)
    summary.to_csv("summary.csv", index=False)
    return summary

# 引数のパース
def parse_args(argv: List[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Minimal splitter for SRA SraRunInfo.csv (flag-style arguments)",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "-i", "--input",
        dest="input_csv", required=True,
        help="Path to SraRunInfo.csv",
    )
    parser.add_argument(
        "-l", "--limit-gb",
        dest="limit_gb", type=float,
        required=True, help="Per-file size limit in GB (1 GB = 1000 MB)",
    )
    parser.add_argument(
        "-p", "--prefix",
        dest="prefix", required=True,
        help="Output filename prefix (e.g., subset)",
    )
    parser.add_argument(
        "-a", "--allow-oversized",
        action="store_true",
        help="Place an oversized Experiment alone in its own file",
    )
    parser.add_argument(
        "-o", "--organism",
        dest="organism", default="human",
        help='Value to write into the "organism" column (default: "human")',
    )
    return parser.parse_args(argv)

# Main
def main(argv: List[str] | None = None) -> int:
    args = parse_args(argv)

    limit_mb = args.limit_gb * 1000.0  # 1 GB = 1000 MB (fixed)

    # Dataframeの読み込み
    df         = pd.read_csv(
        args.input_csv,
        usecols=["Experiment", "Run", "size_MB", "TaxID"],
        dtype={"Experiment": "string", "Run": "string", size_MB: "float", "TaxID": "Int64"},
    )
    # TaxID→organismのフィルタリング
    if args.organism.lower() == "human":
        df = df[df["organism"] == 9606]
    elif args.organism.lower() == "mouse":
        df = df[df["organism"] == 10090]
    elif args.organism.lower() == "rat":
        df = df[df["organism"] == 10116]
    else:
        print(f"Warning: Unrecognized organism '{args.organism}'. No filtering applied.", file=sys.stderr)
    #  TaxID列の削除
    df = df.drop(columns=["TaxID"])

    # Experiment list
    exp_order  = df["Experiment"].drop_duplicates(keep="first").tolist()
    # Experimentの合計サイズ
    exp_totals = df.groupby("Experiment", sort=False)["size_MB"].sum().to_dict()

    # ファイルサイズごとにExperimentの分離
    bins = pack_inorder(exp_order,
                        exp_totals,
                        limit_mb=limit_mb,
                        allow_oversized=args.allow_oversized)

    # Print plan
    print(
        f"Input runs: {df.shape[0]}; Experiments: {len(exp_totals)}; Files: {len(bins)}; "
        f"Limit: {limit_mb:.2f} MB (~{args.limit_gb:.2f} GB) | organism={args.organism}"
    )

    # Show each bin
    for i, b in enumerate(bins, start=1):
        runs = int(df[df["Experiment"].isin(b.experiments)].shape[0])
        print(
            f"  {i:03d}: exps={len(b.experiments):4d} runs={runs:5d} total={b.total_mb/1000.0:.2f} GB ({b.total_mb:.2f} MB)"
        )

    # export
    summary = write_outputs(df, bins, args.prefix, args.organism)
    print("\nSummary written to summary.csv")
    print(summary.to_string(index=False))
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
