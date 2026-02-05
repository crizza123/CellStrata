#!/usr/bin/env python3
"""
download_and_qc.py - Download a specific CELLxGENE Census dataset and perform QC.

This script:
1. Downloads dataset 1b350d0a-4535-4879-beb6-1142f3f94947 from CELLxGENE Census
2. Saves it as .h5ad format
3. Performs QC metrics calculation
4. Generates QC visualizations
5. Applies QC filtering

Usage:
    python download_and_qc.py

Requirements:
    pip install cellxgene-census scanpy matplotlib seaborn

Author: CellStrata Team
Course: BIOINF 576
"""

import os
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

# Create output directories
Path("data/raw").mkdir(parents=True, exist_ok=True)
Path("data/processed").mkdir(parents=True, exist_ok=True)
Path("figures/qc").mkdir(parents=True, exist_ok=True)

print("=" * 60)
print("CellStrata: Dataset Download and QC Pipeline")
print("=" * 60)

# =============================================================================
# Step 1: Download the dataset
# =============================================================================
print("\n[Step 1] Downloading dataset from CELLxGENE Census...")

import cellxgene_census as cxc
import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# Dataset ID to download
DATASET_ID = "1b350d0a-4535-4879-beb6-1142f3f94947"
OUTPUT_RAW = "data/raw/dataset_1b350d0a.h5ad"
OUTPUT_QC = "data/processed/dataset_1b350d0a_qc.h5ad"

# Check if already downloaded
if Path(OUTPUT_RAW).exists():
    print(f"  Dataset already exists at {OUTPUT_RAW}")
    print("  Loading from disk...")
    adata = sc.read_h5ad(OUTPUT_RAW)
else:
    print(f"  Querying Census for dataset_id: {DATASET_ID}")
    
    # Open Census and download
    with cxc.open_soma(census_version="stable") as census:
        # First, let's check the dataset info
        print("  Fetching dataset metadata...")
        
        # Get the experiment
        exp = census["census_data"]["homo_sapiens"]
        
        # Build filter for this specific dataset
        value_filter = f"dataset_id == '{DATASET_ID}'"
        
        print(f"  Filter: {value_filter}")
        print("  Downloading AnnData (this may take a few minutes)...")
        
        # Download as AnnData
        adata = cxc.get_anndata(
            census=census,
            organism="homo_sapiens",
            obs_value_filter=value_filter,
            column_names={
                "obs": [
                    "dataset_id", "donor_id", "assay", "cell_type", 
                    "tissue", "tissue_general", "disease", "sex",
                    "development_stage", "is_primary_data", "suspension_type",
                    "cell_type_ontology_term_id", "tissue_ontology_term_id",
                    "disease_ontology_term_id", "assay_ontology_term_id",
                ],
            },
        )
        
        print(f"  Downloaded {adata.n_obs:,} cells x {adata.n_vars:,} genes")
        
        # Save raw data
        print(f"  Saving to {OUTPUT_RAW}...")
        adata.write_h5ad(OUTPUT_RAW)
        print("  ✓ Raw data saved")

print(f"\n  Dataset shape: {adata.n_obs:,} cells x {adata.n_vars:,} genes")

# Show dataset metadata summary
print("\n  Dataset Metadata Summary:")
print(f"    Tissues: {adata.obs['tissue_general'].unique().tolist()}")
print(f"    Assays: {adata.obs['assay'].unique().tolist()}")
print(f"    Disease: {adata.obs['disease'].unique().tolist()}")
print(f"    Cell types: {adata.obs['cell_type'].nunique()} unique types")

# =============================================================================
# Step 2: Calculate QC Metrics
# =============================================================================
print("\n[Step 2] Calculating QC metrics...")

# Make a copy for QC
adata_qc = adata.copy()

# Store raw counts in a layer
adata_qc.layers['counts'] = adata_qc.X.copy()

# Identify mitochondrial genes
adata_qc.var['mt'] = adata_qc.var_names.str.startswith('MT-')
print(f"  Found {adata_qc.var['mt'].sum()} mitochondrial genes")

# Identify ribosomal genes
adata_qc.var['ribo'] = adata_qc.var_names.str.startswith(('RPS', 'RPL'))
print(f"  Found {adata_qc.var['ribo'].sum()} ribosomal genes")

# Calculate QC metrics
print("  Computing QC metrics...")
sc.pp.calculate_qc_metrics(
    adata_qc, 
    qc_vars=['mt', 'ribo'], 
    percent_top=None, 
    log1p=False, 
    inplace=True
)

# QC metrics are now in adata_qc.obs:
# - n_genes_by_counts: number of genes with positive counts
# - total_counts: total UMI counts per cell
# - pct_counts_mt: percentage of counts from mitochondrial genes
# - pct_counts_ribo: percentage of counts from ribosomal genes

print("\n  QC Metrics Summary:")
print(f"    n_genes_by_counts: median={adata_qc.obs['n_genes_by_counts'].median():.0f}, "
      f"range=[{adata_qc.obs['n_genes_by_counts'].min():.0f}, {adata_qc.obs['n_genes_by_counts'].max():.0f}]")
print(f"    total_counts: median={adata_qc.obs['total_counts'].median():.0f}, "
      f"range=[{adata_qc.obs['total_counts'].min():.0f}, {adata_qc.obs['total_counts'].max():.0f}]")
print(f"    pct_counts_mt: median={adata_qc.obs['pct_counts_mt'].median():.2f}%, "
      f"range=[{adata_qc.obs['pct_counts_mt'].min():.2f}%, {adata_qc.obs['pct_counts_mt'].max():.2f}%]")

# =============================================================================
# Step 3: Generate QC Visualizations
# =============================================================================
print("\n[Step 3] Generating QC visualizations...")

# Set style
sc.settings.set_figure_params(dpi=100, facecolor='white', figsize=(12, 4))

# 3.1: Violin plots of QC metrics
print("  Creating violin plots...")
fig, axes = plt.subplots(1, 4, figsize=(16, 4))

# n_genes violin
sc.pl.violin(adata_qc, 'n_genes_by_counts', ax=axes[0], show=False)
axes[0].set_title('Genes per Cell')
axes[0].set_ylabel('n_genes_by_counts')

# total_counts violin
sc.pl.violin(adata_qc, 'total_counts', ax=axes[1], show=False)
axes[1].set_title('UMI Counts per Cell')
axes[1].set_ylabel('total_counts')

# pct_counts_mt violin
sc.pl.violin(adata_qc, 'pct_counts_mt', ax=axes[2], show=False)
axes[2].set_title('% Mitochondrial')
axes[2].set_ylabel('pct_counts_mt')

# pct_counts_ribo violin
sc.pl.violin(adata_qc, 'pct_counts_ribo', ax=axes[3], show=False)
axes[3].set_title('% Ribosomal')
axes[3].set_ylabel('pct_counts_ribo')

plt.tight_layout()
plt.savefig('figures/qc/qc_violin_plots.png', dpi=150, bbox_inches='tight')
plt.close()
print("    ✓ Saved figures/qc/qc_violin_plots.png")

# 3.2: Scatter plots
print("  Creating scatter plots...")
fig, axes = plt.subplots(1, 3, figsize=(15, 4))

# Total counts vs n_genes
axes[0].scatter(adata_qc.obs['total_counts'], adata_qc.obs['n_genes_by_counts'], 
                alpha=0.1, s=1, c='steelblue')
axes[0].set_xlabel('Total UMI Counts')
axes[0].set_ylabel('Number of Genes')
axes[0].set_title('Counts vs Genes')

# Total counts vs pct_mt
axes[1].scatter(adata_qc.obs['total_counts'], adata_qc.obs['pct_counts_mt'], 
                alpha=0.1, s=1, c='steelblue')
axes[1].set_xlabel('Total UMI Counts')
axes[1].set_ylabel('% Mitochondrial')
axes[1].set_title('Counts vs % Mito')

# n_genes vs pct_mt
axes[2].scatter(adata_qc.obs['n_genes_by_counts'], adata_qc.obs['pct_counts_mt'], 
                alpha=0.1, s=1, c='steelblue')
axes[2].set_xlabel('Number of Genes')
axes[2].set_ylabel('% Mitochondrial')
axes[2].set_title('Genes vs % Mito')

plt.tight_layout()
plt.savefig('figures/qc/qc_scatter_plots.png', dpi=150, bbox_inches='tight')
plt.close()
print("    ✓ Saved figures/qc/qc_scatter_plots.png")

# 3.3: Histograms with thresholds
print("  Creating histograms with threshold lines...")
fig, axes = plt.subplots(2, 2, figsize=(12, 10))

# Define QC thresholds
MIN_GENES = 200
MAX_GENES = 6000
MAX_PCT_MT = 20
MIN_COUNTS = 500

# n_genes histogram
axes[0, 0].hist(adata_qc.obs['n_genes_by_counts'], bins=100, color='steelblue', edgecolor='none')
axes[0, 0].axvline(MIN_GENES, color='red', linestyle='--', label=f'min={MIN_GENES}')
axes[0, 0].axvline(MAX_GENES, color='red', linestyle='--', label=f'max={MAX_GENES}')
axes[0, 0].set_xlabel('Number of Genes')
axes[0, 0].set_ylabel('Number of Cells')
axes[0, 0].set_title('Distribution of Genes per Cell')
axes[0, 0].legend()

# total_counts histogram
axes[0, 1].hist(adata_qc.obs['total_counts'], bins=100, color='steelblue', edgecolor='none')
axes[0, 1].axvline(MIN_COUNTS, color='red', linestyle='--', label=f'min={MIN_COUNTS}')
axes[0, 1].set_xlabel('Total UMI Counts')
axes[0, 1].set_ylabel('Number of Cells')
axes[0, 1].set_title('Distribution of UMI Counts per Cell')
axes[0, 1].legend()

# pct_mt histogram
axes[1, 0].hist(adata_qc.obs['pct_counts_mt'], bins=100, color='steelblue', edgecolor='none')
axes[1, 0].axvline(MAX_PCT_MT, color='red', linestyle='--', label=f'max={MAX_PCT_MT}%')
axes[1, 0].set_xlabel('% Mitochondrial')
axes[1, 0].set_ylabel('Number of Cells')
axes[1, 0].set_title('Distribution of % Mitochondrial')
axes[1, 0].legend()

# Summary stats table
axes[1, 1].axis('off')
stats_text = f"""
QC Thresholds:
─────────────────────────────
Min genes per cell:    {MIN_GENES}
Max genes per cell:    {MAX_GENES}
Max % mitochondrial:   {MAX_PCT_MT}%
Min UMI counts:        {MIN_COUNTS}

Current Dataset:
─────────────────────────────
Total cells:           {adata_qc.n_obs:,}
Total genes:           {adata_qc.n_vars:,}

Cells passing filters:
  n_genes ≥ {MIN_GENES}:        {(adata_qc.obs['n_genes_by_counts'] >= MIN_GENES).sum():,}
  n_genes ≤ {MAX_GENES}:        {(adata_qc.obs['n_genes_by_counts'] <= MAX_GENES).sum():,}
  pct_mt ≤ {MAX_PCT_MT}%:        {(adata_qc.obs['pct_counts_mt'] <= MAX_PCT_MT).sum():,}
  total_counts ≥ {MIN_COUNTS}:   {(adata_qc.obs['total_counts'] >= MIN_COUNTS).sum():,}
"""
axes[1, 1].text(0.1, 0.9, stats_text, transform=axes[1, 1].transAxes, 
                fontsize=11, verticalalignment='top', fontfamily='monospace',
                bbox=dict(boxstyle='round', facecolor='lightgray', alpha=0.5))

plt.tight_layout()
plt.savefig('figures/qc/qc_histograms.png', dpi=150, bbox_inches='tight')
plt.close()
print("    ✓ Saved figures/qc/qc_histograms.png")

# =============================================================================
# Step 4: Apply QC Filtering
# =============================================================================
print("\n[Step 4] Applying QC filters...")

n_cells_before = adata_qc.n_obs
n_genes_before = adata_qc.n_vars

# Filter cells
print(f"  Before filtering: {n_cells_before:,} cells")

# Apply cell filters
sc.pp.filter_cells(adata_qc, min_genes=MIN_GENES)
print(f"    After min_genes={MIN_GENES}: {adata_qc.n_obs:,} cells")

adata_qc = adata_qc[adata_qc.obs['n_genes_by_counts'] <= MAX_GENES, :].copy()
print(f"    After max_genes={MAX_GENES}: {adata_qc.n_obs:,} cells")

adata_qc = adata_qc[adata_qc.obs['pct_counts_mt'] <= MAX_PCT_MT, :].copy()
print(f"    After max_pct_mt={MAX_PCT_MT}%: {adata_qc.n_obs:,} cells")

adata_qc = adata_qc[adata_qc.obs['total_counts'] >= MIN_COUNTS, :].copy()
print(f"    After min_counts={MIN_COUNTS}: {adata_qc.n_obs:,} cells")

# Filter genes (keep genes expressed in at least 3 cells)
MIN_CELLS = 3
sc.pp.filter_genes(adata_qc, min_cells=MIN_CELLS)
print(f"    After min_cells={MIN_CELLS} for genes: {adata_qc.n_vars:,} genes")

n_cells_after = adata_qc.n_obs
n_genes_after = adata_qc.n_vars

print(f"\n  Filtering Summary:")
print(f"    Cells: {n_cells_before:,} → {n_cells_after:,} ({n_cells_before - n_cells_after:,} removed, {100*n_cells_after/n_cells_before:.1f}% retained)")
print(f"    Genes: {n_genes_before:,} → {n_genes_after:,} ({n_genes_before - n_genes_after:,} removed, {100*n_genes_after/n_genes_before:.1f}% retained)")

# =============================================================================
# Step 5: Save QC-filtered data
# =============================================================================
print("\n[Step 5] Saving QC-filtered data...")

# Save
adata_qc.write_h5ad(OUTPUT_QC)
print(f"  ✓ Saved to {OUTPUT_QC}")

# =============================================================================
# Step 6: Generate post-QC summary
# =============================================================================
print("\n[Step 6] Post-QC Summary...")

# Create summary CSV
qc_summary = pd.DataFrame({
    'Metric': [
        'Total cells (raw)',
        'Total genes (raw)',
        'Total cells (after QC)',
        'Total genes (after QC)',
        'Cells removed',
        'Genes removed',
        'Cell retention rate',
        'Gene retention rate',
        'Min genes threshold',
        'Max genes threshold',
        'Max % mito threshold',
        'Min counts threshold',
        'Min cells for gene',
    ],
    'Value': [
        f"{n_cells_before:,}",
        f"{n_genes_before:,}",
        f"{n_cells_after:,}",
        f"{n_genes_after:,}",
        f"{n_cells_before - n_cells_after:,}",
        f"{n_genes_before - n_genes_after:,}",
        f"{100*n_cells_after/n_cells_before:.1f}%",
        f"{100*n_genes_after/n_genes_before:.1f}%",
        str(MIN_GENES),
        str(MAX_GENES),
        f"{MAX_PCT_MT}%",
        str(MIN_COUNTS),
        str(MIN_CELLS),
    ]
})

qc_summary.to_csv('figures/qc/qc_summary.csv', index=False)
print("  ✓ Saved figures/qc/qc_summary.csv")

# Print summary table
print("\n" + "=" * 60)
print("QC COMPLETE")
print("=" * 60)
print(qc_summary.to_string(index=False))
print("=" * 60)

print("\n  Output files:")
print(f"    Raw data:        {OUTPUT_RAW}")
print(f"    QC-filtered:     {OUTPUT_QC}")
print(f"    Violin plots:    figures/qc/qc_violin_plots.png")
print(f"    Scatter plots:   figures/qc/qc_scatter_plots.png")
print(f"    Histograms:      figures/qc/qc_histograms.png")
print(f"    Summary CSV:     figures/qc/qc_summary.csv")

print("\n✓ Pipeline complete!")