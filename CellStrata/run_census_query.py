#!/usr/bin/env python3
"""
Standalone Census query script — works exactly like Untitled-1.ipynb.

Usage:
    python run_census_query.py

Edit the filters below to match your query. Output is a CSV of
dataset_id + cell_count, plus a pandas DataFrame you can inspect.
"""

import cellxgene_census as cxc
import pandas as pd

# ── Filters (edit these) ────────────────────────────────────────────
CENSUS_VERSION = "stable"
ORGANISM = "homo_sapiens"

VALUE_FILTER = (
    "is_primary_data == True"
    " and disease_ontology_term_id == 'PATO:0000461'"  # healthy/normal
    " and tissue_general == 'brain'"
    " and cell_type == 'T cell'"
)

OBS_COLS = [
    "dataset_id",
    "donor_id",
    "assay",
    "sex",
    "tissue",
    "tissue_general",
    "cell_type",
    "disease",
    "development_stage",
]

OUTPUT_CSV = "results.csv"
# ────────────────────────────────────────────────────────────────────

print(f"Connecting to Census ({CENSUS_VERSION})...")

with cxc.open_soma(census_version=CENSUS_VERSION) as census:
    human = census["census_data"][ORGANISM]

    print("Running query...")
    df = (
        human.obs.read(
            value_filter=VALUE_FILTER,
            column_names=OBS_COLS,
        )
        .concat()
        .to_pandas()
    )
    print(f"Got {len(df):,} cells")

    # Dataset summary
    summary = (
        df.groupby("dataset_id", sort=True)
        .size()
        .reset_index(name="cell_count")
    )
    summary.to_csv(OUTPUT_CSV, index=False)
    print(f"\nWrote {OUTPUT_CSV} ({len(summary)} datasets)")
    print(summary.to_string(index=False))

    # Quick stats
    print(f"\nAssays:  {df['assay'].value_counts().to_dict()}")
    print(f"Sex:     {df['sex'].value_counts().to_dict()}")
    print(f"Donors:  {df['donor_id'].nunique()}")

print("\nDone.")
