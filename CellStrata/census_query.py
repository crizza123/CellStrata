#!/usr/bin/env python3
"""
census_query.py - YAML-driven CELLxGENE Census query module for CellStrata.

#Disclaimer AI was used to help write this module. Please review carefully.

This module provides functionality to query and download single-cell RNA-seq
data from the CELLxGENE Census database. It supports:
  - YAML configuration for reproducible queries
  - Label-to-ontology-term-id resolution
  - Multiple output formats: pandas DataFrame, Arrow Table, AnnData, Parquet
    - AnnData for all purpose useage
    - Parquet for metadata summaries 
    - Pandas and Arrow for flexible downstream processing
  - Streaming large queries to avoid OOM errors
  - TileDB context tuning for S3 timeout mitigation

Usage:
    # As a module
    from census_query import load_query_spec_yaml, run_query
    spec = load_query_spec_yaml("config/census_query.yaml")
    result = run_query(spec)

    # As a CLI
    python census_query.py --config config/census_query.yaml

Author: Crizza
Course: BIOINF 576
"""

from __future__ import annotations

import argparse
import logging
import os
from dataclasses import dataclass, field
from pathlib import Path
from typing import (
    Any,
    Dict,
    Iterable,
    Iterator,
    List,
    Literal,
    Optional,
    Sequence,
    Tuple,
    Union,
)

import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
import yaml

import cellxgene_census as cxc

# Optional tiledbsoma for context tuning
try:
    import tiledbsoma as soma
except ImportError:
    soma = None

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
)
logger = logging.getLogger(__name__)


# =============================================================================
# Schema Constants
# =============================================================================
OBS_TERM_ID_COLS: Dict[str, str] = {
    "assay": "assay_ontology_term_id",
    "cell_type": "cell_type_ontology_term_id",
    "development_stage": "development_stage_ontology_term_id",
    "disease": "disease_ontology_term_id",
    "self_reported_ethnicity": "self_reported_ethnicity_ontology_term_id",
    "sex": "sex_ontology_term_id",
    "tissue": "tissue_ontology_term_id",
    "tissue_general": "tissue_general_ontology_term_id",
}

OBS_LABEL_COLS: Dict[str, str] = {
    "assay": "assay",
    "cell_type": "cell_type",
    "development_stage": "development_stage",
    "disease": "disease",
    "self_reported_ethnicity": "self_reported_ethnicity",
    "sex": "sex",
    "tissue": "tissue",
    "tissue_general": "tissue_general",
}


# =============================================================================
# Configuration Dataclasses
# =============================================================================
@dataclass(frozen=True)
class CensusTarget:
    """Target Census dataset specification."""
    census_version: str = "stable"
    organism: str = "homo_sapiens"
    measurement: str = "RNA"


@dataclass(frozen=True)
class ObsFilters:
    """Observation (cell) metadata filters."""
    # Canonical filters
    is_primary_data: Optional[bool] = True
    suspension_type: Optional[str] = "cell"

    # Disease filters (prefer ontology term IDs for stability)
    disease_ontology_term_ids: List[str] = field(default_factory=list)
    disease_labels: List[str] = field(default_factory=list)

    # Sex filters
    sex_labels: List[str] = field(default_factory=list)
    sex_ontology_term_ids: List[str] = field(default_factory=list)

    # Assay filters
    assay_labels: List[str] = field(default_factory=list)
    assay_ontology_term_ids: List[str] = field(default_factory=list)

    # Tissue/organ filters (coarse)
    tissue_general_labels: List[str] = field(default_factory=list)
    tissue_general_ontology_term_ids: List[str] = field(default_factory=list)

    # Tissue filters (fine-grained)
    tissue_labels: List[str] = field(default_factory=list)
    tissue_ontology_term_ids: List[str] = field(default_factory=list)

    # Cell type filters
    cell_type_labels: List[str] = field(default_factory=list)
    cell_type_ontology_term_ids: List[str] = field(default_factory=list)

    # Development stage filters
    development_stage_labels: List[str] = field(default_factory=list)
    development_stage_ontology_term_ids: List[str] = field(default_factory=list)

    # Additional custom filter expression
    extra_value_filter: Optional[str] = None


@dataclass(frozen=True)
class OutputSpec:
    """Output specification."""
    mode: Literal["pandas", "arrow", "anndata", "parquet"] = "pandas"
    outpath: Optional[str] = None
    overwrite: bool = True
    parquet_compression: str = "zstd"


@dataclass(frozen=True)
class QuerySpec:
    """Complete query specification."""
    target: CensusTarget = field(default_factory=CensusTarget)
    obs_filters: ObsFilters = field(default_factory=ObsFilters)

    # Output column configuration
    export_all_non_ontology_obs_columns: bool = False
    obs_columns: List[str] = field(default_factory=list)

    # AnnData-specific options
    var_value_filter: Optional[str] = None
    var_columns: List[str] = field(default_factory=list)

    # TileDB configuration for S3 timeout mitigation
    tiledb_config: Dict[str, Any] = field(default_factory=dict)

    output: OutputSpec = field(default_factory=OutputSpec)


# =============================================================================
# YAML Loading
# =============================================================================
def load_query_spec_yaml(path: Union[str, Path]) -> QuerySpec:
    """
    Load a QuerySpec from a YAML configuration file.

    Args:
        path: Path to the YAML configuration file.

    Returns:
        QuerySpec instance populated from the YAML file.

    Raises:
        FileNotFoundError: If the configuration file doesn't exist.
        yaml.YAMLError: If the YAML is malformed.
    """
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"Configuration file not found: {path}")

    logger.info(f"Loading query specification from: {path}")

    with path.open("r", encoding="utf-8") as f:
        raw = yaml.safe_load(f)

    target = CensusTarget(**raw.get("target", {}))
    obs_filters = ObsFilters(**raw.get("obs_filters", {}))
    output = OutputSpec(**raw.get("output", {}))

    return QuerySpec(
        target=target,
        obs_filters=obs_filters,
        export_all_non_ontology_obs_columns=raw.get(
            "export_all_non_ontology_obs_columns", False
        ),
        obs_columns=raw.get("obs_columns", []) or [],
        var_value_filter=raw.get("var_value_filter"),
        var_columns=raw.get("var_columns", []) or [],
        tiledb_config=raw.get("tiledb_config", {}) or {},
        output=output,
    )


# =============================================================================
# Helper Functions
# =============================================================================
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
    return (
        census["census_info"]["summary_cell_counts"]
        .read(column_names=["category", "label", "ontology_term_id"])
        .concat()
        .to_pandas()
    )


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

    # Development stage filter
    ds_ids = list(obs_filters.development_stage_ontology_term_ids)
    if not ds_ids and obs_filters.development_stage_labels:
        ds_ids = resolve_ontology_term_ids(
            census, "development_stage", obs_filters.development_stage_labels
        )
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


# =============================================================================
# I/O Functions
# =============================================================================
def _ensure_writable_parent(path: Path) -> None:
    """Ensure the parent directory exists and is writable."""
    path.parent.mkdir(parents=True, exist_ok=True)
    if not os.access(path.parent, os.W_OK):
        raise PermissionError(f"Output directory is not writable: {path.parent}")


def _resolve_outpath(outpath: str, overwrite: bool) -> Path:
    """
    Resolve and validate the output path.

    Args:
        outpath: Target output path.
        overwrite: Whether to overwrite existing files.

    Returns:
        Resolved Path object.

    Raises:
        FileExistsError: If file exists and overwrite is False.
        PermissionError: If file cannot be written/overwritten.
    """
    p = Path(outpath)
    _ensure_writable_parent(p)
    if p.exists():
        if overwrite:
            if not os.access(p, os.W_OK):
                raise PermissionError(f"Cannot overwrite existing file: {p}")
            p.unlink()
        else:
            raise FileExistsError(f"Output exists and overwrite=False: {p}")
    return p


def _make_soma_context(tiledb_config: Dict[str, Any]):
    """
    Create a SOMATileDBContext for tuning S3 timeouts/retries.

    Args:
        tiledb_config: TileDB configuration dictionary.

    Returns:
        SOMATileDBContext or None if not needed/available.
    """
    if not tiledb_config:
        return None
    if soma is None:
        raise RuntimeError(
            "tiledbsoma is not importable, but tiledb_config was provided."
        )
    return soma.SOMATileDBContext(tiledb_config=tiledb_config)


def stream_obs_tables(
    exp,
    value_filter: str,
    column_names: Sequence[str],
) -> Iterator[pa.Table]:
    """
    Stream obs results as Arrow Tables.

    Args:
        exp: Census experiment object.
        value_filter: SOMA value_filter string.
        column_names: Columns to include.

    Yields:
        Arrow Table for each batch of results.
    """
    it = exp.obs.read(value_filter=value_filter, column_names=list(column_names))
    for tbl in it:
        if isinstance(tbl, pa.RecordBatch):
            yield pa.Table.from_batches([tbl])
        elif isinstance(tbl, pa.Table):
            yield tbl
        else:
            yield pa.Table.from_batches(tbl.to_batches())


def write_parquet_stream(
    tables: Iterable[pa.Table],
    outpath: Path,
    compression: str = "zstd",
) -> Tuple[int, int]:
    """
    Write an iterable of Arrow Tables to a single Parquet file.

    Args:
        tables: Iterable of Arrow Tables.
        outpath: Output file path.
        compression: Parquet compression codec.

    Returns:
        Tuple of (rows_written, batches_written).
    """
    writer: Optional[pq.ParquetWriter] = None
    rows = 0
    batches = 0

    try:
        for tbl in tables:
            if tbl.num_rows == 0:
                continue
            if writer is None:
                writer = pq.ParquetWriter(str(outpath), tbl.schema, compression=compression)
            writer.write_table(tbl)
            rows += tbl.num_rows
            batches += 1
            if batches % 10 == 0:
                logger.info(f"Progress: {batches} batches, {rows:,} rows written")
    finally:
        if writer is not None:
            writer.close()

    return rows, batches


# =============================================================================
# Main Query Function
# =============================================================================
def run_query(spec: QuerySpec) -> Union[pd.DataFrame, pa.Table, Any, Path]:
    """
    Execute a Census query based on the QuerySpec.

    Args:
        spec: Query specification.

    Returns:
        Depending on output.mode:
        - "pandas": pandas DataFrame
        - "arrow": Arrow Table
        - "anndata": AnnData object
        - "parquet": Path to output Parquet file

    Raises:
        ValueError: If output mode is invalid or required options are missing.
    """
    logger.info(f"Starting Census query: version={spec.target.census_version}, "
                f"organism={spec.target.organism}")

    context = _make_soma_context(spec.tiledb_config)

    with cxc.open_soma(census_version=spec.target.census_version, context=context) as census:
        exp = census["census_data"][spec.target.organism]
        obs_keys = list(exp.obs.keys())
        organ_col = choose_organ_column(obs_keys)

        # Build the value filter
        value_filter = build_obs_value_filter(census, spec.obs_filters)
        logger.info(f"Value filter: {value_filter}")

        # Determine export columns
        export_cols = build_obs_export_columns(
            obs_keys=obs_keys,
            export_all_non_ontology=spec.export_all_non_ontology_obs_columns,
            user_cols=spec.obs_columns,
        )

        # Ensure core metadata columns are present
        for required in ("cell_type", "tissue", organ_col):
            if required in obs_keys and required not in export_cols:
                export_cols.append(required)

        mode = spec.output.mode
        logger.info(f"Output mode: {mode}")

        # Execute based on output mode
        if mode == "anndata":
            column_names = {"obs": export_cols}
            if spec.var_columns:
                column_names["var"] = spec.var_columns

            organism_name = (
                "Homo sapiens"
                if spec.target.organism == "homo_sapiens"
                else "Mus musculus"
            )

            adata = cxc.get_anndata(
                census=census,
                organism=organism_name,
                obs_value_filter=value_filter or None,
                var_value_filter=spec.var_value_filter or None,
                column_names=column_names,
            )
            logger.info(f"Retrieved AnnData: {adata.n_obs} cells x {adata.n_vars} genes")
            return adata

        if mode == "pandas":
            df = (
                exp.obs.read(value_filter=value_filter or None, column_names=export_cols)
                .concat()
                .to_pandas()
            )
            logger.info(f"Retrieved DataFrame: {len(df):,} rows")
            return df

        if mode == "arrow":
            tables = list(
                stream_obs_tables(exp, value_filter=value_filter, column_names=export_cols)
            )
            if not tables:
                return pa.table({c: [] for c in export_cols})
            result = pa.concat_tables(tables, promote=True)
            logger.info(f"Retrieved Arrow Table: {result.num_rows:,} rows")
            return result

        if mode == "parquet":
            if not spec.output.outpath:
                raise ValueError("output.outpath is required for mode='parquet'")
            outpath = _resolve_outpath(spec.output.outpath, spec.output.overwrite)

            tables = stream_obs_tables(
                exp, value_filter=value_filter, column_names=export_cols
            )
            rows, batches = write_parquet_stream(
                tables, outpath, compression=spec.output.parquet_compression
            )
            logger.info(f"Wrote Parquet: {rows:,} rows in {batches} batches to {outpath}")
            return outpath

        raise ValueError(f"Unknown output mode: {mode}")


# =============================================================================
# CLI Entry Point
# =============================================================================
def main():
    """Command-line interface for census_query."""
    parser = argparse.ArgumentParser(
        description="Query CELLxGENE Census using YAML configuration.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Example:
    python census_query.py --config census_query.yaml

The YAML configuration file should specify:
    - target: Census version, organism, measurement
    - obs_filters: Cell filtering criteria
    - output: Output mode and path
        """,
    )
    parser.add_argument(
        "--config",
        "-c",
        required=True,
        help="Path to YAML configuration file",
    )
    parser.add_argument(
        "--verbose",
        "-v",
        action="store_true",
        help="Enable verbose logging",
    )

    args = parser.parse_args()

    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)

    try:
        spec = load_query_spec_yaml(args.config)
        result = run_query(spec)

        if isinstance(result, pd.DataFrame):
            print(f"\nQuery returned DataFrame with {len(result):,} rows")
            print(f"Columns: {list(result.columns)}")
            print(f"\nFirst 5 rows:\n{result.head()}")
        elif isinstance(result, pa.Table):
            print(f"\nQuery returned Arrow Table with {result.num_rows:,} rows")
        elif isinstance(result, Path):
            print(f"\nQuery output written to: {result}")
        else:
            print(f"\nQuery returned: {type(result)}")

    except Exception as e:
        logger.error(f"Query failed: {e}")
        raise SystemExit(1)


if __name__ == "__main__":
    main()