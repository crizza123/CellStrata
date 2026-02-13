#!/usr/bin/env python3
"""Run a Census query from a YAML config file."""
from __future__ import annotations

import sys
from pathlib import Path

# When running standalone (not installed via pip), ensure the repo root is on
# sys.path so that "from CellStrata.census..." resolves.
_repo_root = str(Path(__file__).resolve().parent.parent.parent)
if _repo_root not in sys.path:
    sys.path.insert(0, _repo_root)

import argparse

from CellStrata.census.config import load_config
from CellStrata.census.query import fetch_obs_df, summarize_by_dataset
from CellStrata.census.outputs import write_summary_csv, write_celllevel_csv


def main():
    ap = argparse.ArgumentParser(description="Run a CellStrata Census query")
    ap.add_argument("--config", required=True, help="Path to YAML config file")
    args = ap.parse_args()

    cfg = load_config(args.config)

    df, value_filter = fetch_obs_df(cfg)
    print("VALUE_FILTER:", value_filter)
    print(f"Got {len(df):,} cells")

    summary = summarize_by_dataset(df)
    write_summary_csv(summary, cfg.output.summary_csv)
    print(f"Wrote {cfg.output.summary_csv} ({len(summary)} datasets)")

    if cfg.output.write_celllevel and cfg.output.celllevel_csv:
        write_celllevel_csv(df, cfg.output.celllevel_csv)
        print(f"Wrote {cfg.output.celllevel_csv} ({len(df):,} cells)")


if __name__ == "__main__":
    main()
