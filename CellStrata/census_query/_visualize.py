"""
_visualize.py - Metadata visualization for filtered Census query results.

This module provides plotting functions for exploring cell metadata
returned by run_query() in pandas mode. Each function accepts a DataFrame
and returns a matplotlib Axes so plots can be further customized.

Requires: matplotlib, seaborn
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import List, Optional, Sequence, Union

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

logger = logging.getLogger(__name__)

# Consistent palette across all plots
_PALETTE = "Set2"


def _validate_column(df: pd.DataFrame, col: str) -> None:
    """Raise KeyError if *col* is missing from *df*."""
    if col not in df.columns:
        raise KeyError(
            f"Column '{col}' not found in DataFrame. "
            f"Available columns: {list(df.columns)}"
        )


def _top_n(series: pd.Series, n: int) -> pd.Series:
    """Return value_counts limited to the top *n* entries."""
    return series.value_counts().head(n)


# ── Individual plot functions ─────────────────────────────────────────


def plot_cell_type_counts(
    df: pd.DataFrame,
    *,
    top_n: int = 20,
    ax: Optional[plt.Axes] = None,
) -> plt.Axes:
    """Horizontal bar chart of cell-type frequencies.

    Args:
        df: Metadata DataFrame (must contain ``cell_type`` column).
        top_n: Number of most-frequent cell types to show.
        ax: Optional existing Axes.

    Returns:
        The matplotlib Axes with the plot.
    """
    _validate_column(df, "cell_type")
    counts = _top_n(df["cell_type"], top_n).sort_values()

    if ax is None:
        _, ax = plt.subplots(figsize=(8, max(4, len(counts) * 0.35)))

    counts.plot.barh(ax=ax, color=sns.color_palette(_PALETTE, len(counts)))
    ax.set_xlabel("Number of cells")
    ax.set_title(f"Top {min(top_n, len(counts))} cell types")
    ax.ticklabel_format(axis="x", style="sci", scilimits=(0, 3))
    return ax


def plot_tissue_counts(
    df: pd.DataFrame,
    *,
    column: str = "tissue_general",
    top_n: int = 20,
    ax: Optional[plt.Axes] = None,
) -> plt.Axes:
    """Horizontal bar chart of tissue / organ frequencies.

    Args:
        df: Metadata DataFrame.
        column: Column name to plot (``tissue_general`` or ``tissue``).
        top_n: Number of entries to show.
        ax: Optional existing Axes.

    Returns:
        The matplotlib Axes.
    """
    _validate_column(df, column)
    counts = _top_n(df[column], top_n).sort_values()

    if ax is None:
        _, ax = plt.subplots(figsize=(8, max(4, len(counts) * 0.35)))

    counts.plot.barh(ax=ax, color=sns.color_palette(_PALETTE, len(counts)))
    ax.set_xlabel("Number of cells")
    ax.set_title(f"Cells per {column.replace('_', ' ')} (top {min(top_n, len(counts))})")
    ax.ticklabel_format(axis="x", style="sci", scilimits=(0, 3))
    return ax


def plot_sex_distribution(
    df: pd.DataFrame,
    *,
    ax: Optional[plt.Axes] = None,
) -> plt.Axes:
    """Bar chart of sex distribution.

    Args:
        df: Metadata DataFrame (must contain ``sex`` column).
        ax: Optional existing Axes.

    Returns:
        The matplotlib Axes.
    """
    _validate_column(df, "sex")
    counts = df["sex"].value_counts()

    if ax is None:
        _, ax = plt.subplots(figsize=(5, 4))

    counts.plot.bar(ax=ax, color=sns.color_palette(_PALETTE, len(counts)))
    ax.set_ylabel("Number of cells")
    ax.set_title("Sex distribution")
    ax.tick_params(axis="x", rotation=0)
    return ax


def plot_assay_counts(
    df: pd.DataFrame,
    *,
    ax: Optional[plt.Axes] = None,
) -> plt.Axes:
    """Bar chart of sequencing assay frequencies.

    Args:
        df: Metadata DataFrame (must contain ``assay`` column).
        ax: Optional existing Axes.

    Returns:
        The matplotlib Axes.
    """
    _validate_column(df, "assay")
    counts = df["assay"].value_counts()

    if ax is None:
        _, ax = plt.subplots(figsize=(6, 4))

    counts.plot.bar(ax=ax, color=sns.color_palette(_PALETTE, len(counts)))
    ax.set_ylabel("Number of cells")
    ax.set_title("Assay distribution")
    ax.tick_params(axis="x", rotation=30)
    plt.setp(ax.get_xticklabels(), ha="right")
    return ax


def plot_disease_counts(
    df: pd.DataFrame,
    *,
    top_n: int = 15,
    ax: Optional[plt.Axes] = None,
) -> plt.Axes:
    """Horizontal bar chart of disease annotation frequencies.

    Args:
        df: Metadata DataFrame (must contain ``disease`` column).
        top_n: Number of entries to show.
        ax: Optional existing Axes.

    Returns:
        The matplotlib Axes.
    """
    _validate_column(df, "disease")
    counts = _top_n(df["disease"], top_n).sort_values()

    if ax is None:
        _, ax = plt.subplots(figsize=(7, max(3, len(counts) * 0.35)))

    counts.plot.barh(ax=ax, color=sns.color_palette(_PALETTE, len(counts)))
    ax.set_xlabel("Number of cells")
    ax.set_title(f"Disease annotations (top {min(top_n, len(counts))})")
    return ax


def plot_dataset_contribution(
    df: pd.DataFrame,
    *,
    top_n: int = 15,
    ax: Optional[plt.Axes] = None,
) -> plt.Axes:
    """Horizontal bar chart showing how many cells come from each dataset.

    Args:
        df: Metadata DataFrame (must contain ``dataset_id`` column).
        top_n: Number of datasets to show.
        ax: Optional existing Axes.

    Returns:
        The matplotlib Axes.
    """
    _validate_column(df, "dataset_id")
    counts = _top_n(df["dataset_id"], top_n).sort_values()

    if ax is None:
        _, ax = plt.subplots(figsize=(8, max(4, len(counts) * 0.35)))

    counts.plot.barh(ax=ax, color=sns.color_palette(_PALETTE, len(counts)))
    ax.set_xlabel("Number of cells")
    ax.set_title(f"Cells per dataset (top {min(top_n, len(counts))})")
    return ax


def plot_development_stage_counts(
    df: pd.DataFrame,
    *,
    top_n: int = 15,
    ax: Optional[plt.Axes] = None,
) -> plt.Axes:
    """Bar chart of development stage frequencies.

    Args:
        df: Metadata DataFrame (must contain ``development_stage`` column).
        top_n: Number of entries to show.
        ax: Optional existing Axes.

    Returns:
        The matplotlib Axes.
    """
    _validate_column(df, "development_stage")
    counts = _top_n(df["development_stage"], top_n)

    if ax is None:
        _, ax = plt.subplots(figsize=(7, 4))

    counts.plot.bar(ax=ax, color=sns.color_palette(_PALETTE, len(counts)))
    ax.set_ylabel("Number of cells")
    ax.set_title(f"Development stages (top {min(top_n, len(counts))})")
    ax.tick_params(axis="x", rotation=40)
    plt.setp(ax.get_xticklabels(), ha="right")
    return ax


def plot_donors_per_dataset(
    df: pd.DataFrame,
    *,
    top_n: int = 15,
    ax: Optional[plt.Axes] = None,
) -> plt.Axes:
    """Bar chart of unique donor count per dataset.

    Args:
        df: Metadata DataFrame (must contain ``dataset_id`` and ``donor_id``).
        top_n: Number of datasets to show.
        ax: Optional existing Axes.

    Returns:
        The matplotlib Axes.
    """
    _validate_column(df, "dataset_id")
    _validate_column(df, "donor_id")

    donor_counts = (
        df.groupby("dataset_id")["donor_id"]
        .nunique()
        .sort_values(ascending=False)
        .head(top_n)
    )

    if ax is None:
        _, ax = plt.subplots(figsize=(8, max(4, len(donor_counts) * 0.35)))

    donor_counts.sort_values().plot.barh(
        ax=ax, color=sns.color_palette(_PALETTE, len(donor_counts))
    )
    ax.set_xlabel("Unique donors")
    ax.set_title(f"Donors per dataset (top {min(top_n, len(donor_counts))})")
    return ax


def plot_cell_type_by_tissue(
    df: pd.DataFrame,
    *,
    top_cell_types: int = 10,
    top_tissues: int = 8,
    tissue_column: str = "tissue_general",
    ax: Optional[plt.Axes] = None,
) -> plt.Axes:
    """Heatmap of cell-type counts across tissues.

    Rows = cell types, columns = tissues. Colour intensity encodes the
    log10(count + 1) so rare populations are still visible.

    Args:
        df: Metadata DataFrame.
        top_cell_types: Number of most-frequent cell types to include.
        top_tissues: Number of most-frequent tissues to include.
        tissue_column: Column name for tissue grouping.
        ax: Optional existing Axes.

    Returns:
        The matplotlib Axes.
    """
    _validate_column(df, "cell_type")
    _validate_column(df, tissue_column)

    top_ct = df["cell_type"].value_counts().head(top_cell_types).index
    top_tis = df[tissue_column].value_counts().head(top_tissues).index

    sub = df[df["cell_type"].isin(top_ct) & df[tissue_column].isin(top_tis)]
    ct = pd.crosstab(sub["cell_type"], sub[tissue_column])
    ct = ct.reindex(index=top_ct, columns=top_tis, fill_value=0)
    ct_log = np.log10(ct + 1)

    if ax is None:
        _, ax = plt.subplots(
            figsize=(max(6, top_tissues * 0.9), max(4, top_cell_types * 0.45))
        )

    sns.heatmap(
        ct_log,
        annot=ct.values,
        fmt="d",
        cmap="YlOrRd",
        linewidths=0.5,
        ax=ax,
        cbar_kws={"label": "log10(count + 1)"},
    )
    ax.set_title("Cell types across tissues")
    ax.set_ylabel("")
    ax.set_xlabel("")
    return ax


# ── Summary dashboard ─────────────────────────────────────────────────


def plot_metadata_summary(
    df: pd.DataFrame,
    *,
    top_n: int = 15,
    tissue_column: str = "tissue_general",
    save_path: Optional[Union[str, Path]] = None,
    dpi: int = 150,
) -> plt.Figure:
    """Multi-panel dashboard summarising the key metadata distributions.

    Creates a 3x2 grid with:
        1. Cell type counts
        2. Tissue counts
        3. Sex distribution
        4. Assay distribution
        5. Dataset contribution
        6. Cell-type-by-tissue heatmap

    Args:
        df: Metadata DataFrame returned by ``run_query()`` in pandas mode.
        top_n: Number of top entries for bar charts.
        tissue_column: Column name for tissue grouping.
        save_path: If provided, saves the figure to this path.
        dpi: Resolution for the saved figure.

    Returns:
        The matplotlib Figure.
    """
    fig, axes = plt.subplots(3, 2, figsize=(16, 18))
    fig.suptitle(
        f"Metadata summary  —  {len(df):,} cells, "
        f"{df['donor_id'].nunique():,} donors, "
        f"{df['dataset_id'].nunique():,} datasets"
        if {"donor_id", "dataset_id"}.issubset(df.columns)
        else f"Metadata summary  —  {len(df):,} cells",
        fontsize=14,
        fontweight="bold",
        y=0.98,
    )

    # 1. Cell type counts
    if "cell_type" in df.columns:
        plot_cell_type_counts(df, top_n=top_n, ax=axes[0, 0])
    else:
        axes[0, 0].set_visible(False)

    # 2. Tissue counts
    if tissue_column in df.columns:
        plot_tissue_counts(df, column=tissue_column, top_n=top_n, ax=axes[0, 1])
    else:
        axes[0, 1].set_visible(False)

    # 3. Sex distribution
    if "sex" in df.columns:
        plot_sex_distribution(df, ax=axes[1, 0])
    else:
        axes[1, 0].set_visible(False)

    # 4. Assay distribution
    if "assay" in df.columns:
        plot_assay_counts(df, ax=axes[1, 1])
    else:
        axes[1, 1].set_visible(False)

    # 5. Dataset contribution
    if "dataset_id" in df.columns:
        plot_dataset_contribution(df, top_n=top_n, ax=axes[2, 0])
    else:
        axes[2, 0].set_visible(False)

    # 6. Cell-type-by-tissue heatmap
    if "cell_type" in df.columns and tissue_column in df.columns:
        plot_cell_type_by_tissue(
            df,
            top_cell_types=min(top_n, 10),
            top_tissues=8,
            tissue_column=tissue_column,
            ax=axes[2, 1],
        )
    else:
        axes[2, 1].set_visible(False)

    fig.tight_layout(rect=[0, 0, 1, 0.96])

    if save_path is not None:
        fig.savefig(save_path, dpi=dpi, bbox_inches="tight")
        logger.info(f"Saved summary figure to {save_path}")

    return fig
