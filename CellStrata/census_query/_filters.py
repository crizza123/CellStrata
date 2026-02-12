"""
_filters.py - Filter building and ontology term resolution.

This module can be used independently of query execution. Import it directly
when you need to build SOMA filters or resolve ontology labels without
running a full query.
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import (
    List,
    Optional,
    Sequence,
)

import pandas as pd

from ._schema import OBS_TERM_ID_COLS, ObsFilters

logger = logging.getLogger(__name__)


def _format_in_list(values: Sequence[str]) -> str:
    """
    Format a list of values for SOMA value_filter syntax.

    Args:
        values: List of string values.

    Returns:
        SOMA-compatible list literal: ['a', 'b', 'c']
    """
    escaped = [v.replace("\\", "\\\\").replace("'", "\\'") for v in values]
    return "[" + ", ".join([f"'{v}'" for v in escaped]) + "]"


def _get_summary_cell_counts(census) -> pd.DataFrame:
    """
    Load the summary_cell_counts table from Census.

    Args:
        census: Open Census connection.

    Returns:
        DataFrame with category, label, and ontology_term_id columns.
    """
    logger.info("Reading summary_cell_counts from Census...")
    result = (
        census["census_info"]["summary_cell_counts"]
        .read(column_names=["category", "label", "ontology_term_id"])
        .concat()
        .to_pandas()
    )
    logger.info(f"Loaded {len(result):,} summary rows")
    return result


def resolve_ontology_term_ids(
    census,
    category: str,
    labels: Sequence[str],
) -> List[str]:
    """
    Resolve human-readable labels to ontology term IDs.

    Args:
        census: Open Census connection.
        category: The category to look up (e.g., "disease", "assay").
        labels: List of labels to resolve.

    Returns:
        List of ontology term IDs in the same order as input labels.

    Raises:
        ValueError: If any label cannot be resolved.
    """
    if not labels:
        return []

    scc = _get_summary_cell_counts(census)
    subset = scc[(scc["category"] == category) & (scc["label"].isin(list(labels)))]
    found = dict(zip(subset["label"], subset["ontology_term_id"]))

    missing = [x for x in labels if x not in found]
    if missing:
        raise ValueError(
            f"Could not resolve ontology_term_id for {category} label(s): {missing}"
        )

    return [found[x] for x in labels]


def choose_organ_column(obs_keys: Sequence[str]) -> str:
    """
    Determine the best column name for organ/tissue filtering.

    Prefers 'tissue_general' (coarse organ) over 'tissue' (fine-grained).

    Args:
        obs_keys: Available column names in obs schema.

    Returns:
        The column name to use for organ filtering.

    Raises:
        ValueError: If neither tissue_general nor tissue is available.
    """
    if "tissue_general" in obs_keys:
        return "tissue_general"
    if "tissue" in obs_keys:
        return "tissue"
    raise ValueError("Neither tissue_general nor tissue found in obs schema keys.")


def detect_inventory_dir(user_hint: Optional[str] = None) -> Path:
    """
    Locate the directory containing development stage inventory CSVs.

    Looks for a folder containing:
      - year_old_stage_pairs_counts.csv
      - development_stage_pairs_counts.csv

    Args:
        user_hint: Optional explicit path (string) to the directory.

    Returns:
        Path to the inventory directory.

    Raises:
        FileNotFoundError: If the directory cannot be found.
    """
    required = ["year_old_stage_pairs_counts.csv", "development_stage_pairs_counts.csv"]

    candidates: List[Path] = []
    if user_hint is not None:
        candidates.append(Path(user_hint))

    candidates.extend([
        Path.home() / "census_exports" / "dev_stage_inventory",
        Path("/mnt/data"),
    ])

    for d in candidates:
        if d.is_dir() and all((d / f).exists() for f in required):
            return d

    raise FileNotFoundError(
        "Could not find inventory CSVs. Set inventory_dir in obs_filters to the "
        "folder containing:\n"
        "  - year_old_stage_pairs_counts.csv\n"
        "  - development_stage_pairs_counts.csv"
    )


def build_stage_ids_ge15(inventory_dir: Path, min_age_years: int = 15) -> List[str]:
    """
    Build a development_stage_ontology_term_id allow-list for age >= min_age_years
    using local CSV inventories.

    Includes all N-year-old stages where N >= min_age_years, plus standard
    adult stage labels (decades, age bins, etc.).

    Args:
        inventory_dir: Path to directory with inventory CSVs.
        min_age_years: Minimum age in years.

    Returns:
        Sorted list of development_stage_ontology_term_id strings.
    """
    year_pairs = pd.read_csv(inventory_dir / "year_old_stage_pairs_counts.csv")
    pairs = pd.read_csv(inventory_dir / "development_stage_pairs_counts.csv")

    ids = set(
        year_pairs.loc[
            year_pairs["age_years"] >= min_age_years,
            "development_stage_ontology_term_id",
        ]
        .astype(str)
        .tolist()
    )

    extra_labels = [
        "15-19 year-old",
        "adult stage",
        "young adult stage",
        "prime adult stage",
        "late adult stage",
        "middle aged stage",
        "third decade stage",
        "fourth decade stage",
        "fifth decade stage",
        "sixth decade stage",
        "seventh decade stage",
        "eighth decade stage",
        "ninth decade stage",
        "60-79 year-old stage",
        "80 year-old and over stage",
        "90 year-old and over stage",
    ]

    extra_ids = (
        pairs.loc[
            pairs["development_stage"].isin(extra_labels),
            "development_stage_ontology_term_id",
        ]
        .astype(str)
        .tolist()
    )
    ids.update(extra_ids)

    ids.discard("unknown")
    logger.info(
        f"Inventory age filter: {len(ids)} development stage IDs for age >= {min_age_years}"
    )
    return sorted(ids)


def build_obs_value_filter(census, obs_filters: ObsFilters) -> str:
    """
    Build a SOMA value_filter string from ObsFilters.

    Automatically resolves labels to ontology term IDs for stable queries.

    Args:
        census: Open Census connection (for label resolution).
        obs_filters: Filter specification.

    Returns:
        SOMA-compatible value_filter string.
    """
    parts: List[str] = []

    # Primary data and suspension type
    if obs_filters.is_primary_data is not None:
        parts.append(f"is_primary_data == {bool(obs_filters.is_primary_data)}")
    if obs_filters.suspension_type:
        parts.append(f"suspension_type == '{obs_filters.suspension_type}'")

    # Disease filter
    disease_ids = list(obs_filters.disease_ontology_term_ids)
    if not disease_ids and obs_filters.disease_labels:
        disease_ids = resolve_ontology_term_ids(
            census, "disease", obs_filters.disease_labels
        )
    if disease_ids:
        parts.append(f"{OBS_TERM_ID_COLS['disease']} in {_format_in_list(disease_ids)}")

    # Sex filter
    sex_ids = list(obs_filters.sex_ontology_term_ids)
    if not sex_ids and obs_filters.sex_labels:
        sex_ids = resolve_ontology_term_ids(census, "sex", obs_filters.sex_labels)
    if sex_ids:
        parts.append(f"{OBS_TERM_ID_COLS['sex']} in {_format_in_list(sex_ids)}")

    # Assay filter
    assay_ids = list(obs_filters.assay_ontology_term_ids)
    if not assay_ids and obs_filters.assay_labels:
        assay_ids = resolve_ontology_term_ids(census, "assay", obs_filters.assay_labels)
    if assay_ids:
        parts.append(f"{OBS_TERM_ID_COLS['assay']} in {_format_in_list(assay_ids)}")

    # Tissue general (coarse organ) filter
    tg_ids = list(obs_filters.tissue_general_ontology_term_ids)
    if not tg_ids and obs_filters.tissue_general_labels:
        tg_ids = resolve_ontology_term_ids(
            census, "tissue_general", obs_filters.tissue_general_labels
        )
    if tg_ids:
        parts.append(
            f"{OBS_TERM_ID_COLS['tissue_general']} in {_format_in_list(tg_ids)}"
        )

    # Tissue (fine-grained) filter
    tissue_ids = list(obs_filters.tissue_ontology_term_ids)
    if not tissue_ids and obs_filters.tissue_labels:
        tissue_ids = resolve_ontology_term_ids(
            census, "tissue", obs_filters.tissue_labels
        )
    if tissue_ids:
        parts.append(f"{OBS_TERM_ID_COLS['tissue']} in {_format_in_list(tissue_ids)}")

    # Cell type filter
    ct_ids = list(obs_filters.cell_type_ontology_term_ids)
    if not ct_ids and obs_filters.cell_type_labels:
        ct_ids = resolve_ontology_term_ids(
            census, "cell_type", obs_filters.cell_type_labels
        )
    if ct_ids:
        parts.append(f"{OBS_TERM_ID_COLS['cell_type']} in {_format_in_list(ct_ids)}")

    # Development stage filter — three sources (in priority order):
    # 1. Explicit ontology term IDs
    # 2. Label resolution
    # 3. Inventory-based age filtering (CSV directory)
    ds_ids = list(obs_filters.development_stage_ontology_term_ids)
    if not ds_ids and obs_filters.development_stage_labels:
        ds_ids = resolve_ontology_term_ids(
            census, "development_stage", obs_filters.development_stage_labels
        )
    if not ds_ids and obs_filters.inventory_dir:
        inv_dir = detect_inventory_dir(obs_filters.inventory_dir)
        ds_ids = build_stage_ids_ge15(inv_dir, obs_filters.min_age_years)
    if ds_ids:
        parts.append(
            f"{OBS_TERM_ID_COLS['development_stage']} in {_format_in_list(ds_ids)}"
        )

    # Extra custom filter
    if obs_filters.extra_value_filter:
        parts.append(f"({obs_filters.extra_value_filter})")

    return " and ".join(parts) if parts else ""


def _default_obs_columns() -> List[str]:
    """Return the default set of obs columns to export."""
    return [
        "dataset_id",
        "donor_id",
        "assay",
        "cell_type",
        "tissue",
        "tissue_general",
        "disease",
        "development_stage",
        "sex",
        "self_reported_ethnicity",
        "is_primary_data",
        "suspension_type",
    ]


def build_obs_export_columns(
    obs_keys: Sequence[str],
    export_all_non_ontology: bool,
    user_cols: Sequence[str],
) -> List[str]:
    """
    Determine which obs columns to include in the output.

    Args:
        obs_keys: Available columns in the obs schema.
        export_all_non_ontology: If True, export all non-ontology columns.
        user_cols: User-specified columns to export.

    Returns:
        List of column names to export.

    Raises:
        ValueError: If requested columns are not in the schema.
    """
    keys = list(obs_keys)

    if export_all_non_ontology:
        return [
            c for c in keys if (c != "soma_joinid" and "ontology_term_id" not in c)
        ]

    cols = list(user_cols) if user_cols else _default_obs_columns()
    missing = [c for c in cols if c not in keys]
    if missing:
        raise ValueError(f"Requested obs columns not present in schema: {missing}")
    return cols
