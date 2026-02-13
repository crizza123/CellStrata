#!/usr/bin/env python3
from __future__ import annotations
import argparse

from cellstrata.census.config import load_config
from cellstrata.census.query import fetch_obs_df, summarize_by_dataset
from cellstrata.census.outputs import write_summary_csv, write_celllevel_csv


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True)
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
