#!/usr/bin/env python3
"""Summarize imputation F-scores from result files."""

import argparse
import glob
import os

import pandas as pd


def parse_args():
    parser = argparse.ArgumentParser(description="Summarize imputation F-scores.")
    parser.add_argument(
        "--root",
        default=".",
        help="Root directory to search recursively for result files.",
    )
    parser.add_argument(
        "--pattern",
        default="**/imputed_per_sample_results.txt",
        help="Glob pattern under --root.",
    )
    parser.add_argument(
        "--path-token",
        default="bin_1",
        help="Only include files whose parent path contains this token.",
    )
    parser.add_argument(
        "--score-column",
        default="F-score",
        help="Column containing per-sample scores.",
    )
    return parser.parse_args()


def main():
    args = parse_args()
    pattern = os.path.join(args.root, args.pattern)
    file_paths = glob.glob(pattern, recursive=True)
    scores = []

    for file_path in file_paths:
        if args.path_token and args.path_token not in file_path:
            continue
        try:
            df = pd.read_csv(file_path, sep=r"\s+", engine="python")
        except (pd.errors.EmptyDataError, OSError) as exc:
            print(f"Skipping file {file_path}: {exc}")
            continue
        if args.score_column in df.columns:
            scores.extend(pd.to_numeric(df[args.score_column], errors="coerce").dropna().tolist())

    if scores:
        average = sum(scores) / len(scores)
        print(f"Average {args.score_column} across {len(scores)} values: {average:.4f}")
    else:
        print(f"No {args.score_column} values found.")


if __name__ == "__main__":
    main()
