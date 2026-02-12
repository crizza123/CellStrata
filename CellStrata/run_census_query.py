#!/usr/bin/env python3
"""
Standalone Census query script — works exactly like Untitled-1.ipynb.

Usage:
    python run_census_query.py

Edit the filters below to match your query. Output is a CSV of
dataset_id + cell_count, plus a pandas DataFrame you can inspect.
"""

import cellxgene_census as cxc
import matplotlib.pyplot as plt
import pandas as pd

from CellStrata.census_query._visualize import (
    plot_assay_counts,
    plot_dataset_contribution,
    plot_donors_per_dataset,
    plot_sex_distribution,
)

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

    # ── Visualization ─────────────────────────────────────────────
    SAVE_PATH = "query_summary.png"

    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    fig.suptitle(
        f"Census query — {len(df):,} cells, "
        f"{df['donor_id'].nunique():,} donors, "
        f"{df['dataset_id'].nunique():,} datasets",
        fontsize=14,
        fontweight="bold",
    )

    plot_assay_counts(df, ax=axes[0, 0])
    plot_sex_distribution(df, ax=axes[0, 1])
    plot_dataset_contribution(df, ax=axes[1, 0])
    plot_donors_per_dataset(df, ax=axes[1, 1])

    fig.tight_layout(rect=[0, 0, 1, 0.95])
    fig.savefig(SAVE_PATH, dpi=150, bbox_inches="tight")
    print(f"\nSaved visualization to {SAVE_PATH}")

print("\nDone.")
