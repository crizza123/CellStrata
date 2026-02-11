"""
_runner.py - Query orchestration and CLI entry point.

This module ties together _schema, _filters, and _visualize to execute
Census queries and optionally produce visualizations.
"""

from __future__ import annotations

import argparse
import logging
import os
from pathlib import Path
from typing import (
    Any,
    Dict,
    List,
    Optional,
    Union,
)

import pandas as pd

import cellxgene_census as cxc

from ._schema import QuerySpec, load_query_spec_yaml
from ._filters import (
    build_obs_value_filter,
    build_obs_export_columns,
    choose_organ_column,
)

# Optional tiledbsoma for context tuning
try:
    import tiledbsoma as soma
except ImportError:
    soma = None

logger = logging.getLogger(__name__)


# ── Helpers (inlined from former _io.py) ──────────────────────────────


def _make_soma_context(tiledb_config: Dict[str, Any]):
    """Create a SOMATileDBContext for tuning S3 timeouts/retries."""
    if not tiledb_config:
        return None
    if soma is None:
        raise RuntimeError(
            "tiledbsoma is not importable, but tiledb_config was provided."
        )
    return soma.SOMATileDBContext(tiledb_config=tiledb_config)


def _resolve_outpath(outpath: str, overwrite: bool) -> Path:
    """Resolve and validate an output file path."""
    p = Path(outpath)
    p.parent.mkdir(parents=True, exist_ok=True)
    if not os.access(p.parent, os.W_OK):
        raise PermissionError(f"Output directory is not writable: {p.parent}")
    if p.exists():
        if overwrite:
            if not os.access(p, os.W_OK):
                raise PermissionError(f"Cannot overwrite existing file: {p}")
            p.unlink()
        else:
            raise FileExistsError(f"Output exists and overwrite=False: {p}")
    return p


# ── Query execution ──────────────────────────────────────────────────


def run_query(spec: QuerySpec) -> Union[pd.DataFrame, Any, List[str]]:
    """
    Execute a Census query based on the QuerySpec.

    Args:
        spec: Query specification.

    Returns:
        Depending on output.mode:
        - "pandas": pandas DataFrame
        - "anndata": AnnData object
        - "dataset_list": sorted list of unique dataset ID strings

    Raises:
        ValueError: If output mode is invalid or required options are missing.
    """
    logger.info(f"Starting Census query: version={spec.target.census_version}, "
                f"organism={spec.target.organism}")

    context = _make_soma_context(spec.tiledb_config)

    with cxc.open_soma(census_version=spec.target.census_version, context=context) as census:
        exp = census["census_data"][spec.target.organism]
        obs_keys = list(exp.obs.keys())

        # Build the value filter
        value_filter = build_obs_value_filter(census, spec.obs_filters)
        logger.info(f"Value filter: {value_filter}")

        mode = spec.output.mode
        logger.info(f"Output mode: {mode}")

        # dataset_list mode only needs dataset_id — skip column/organ setup
        if mode == "dataset_list":
            df = (
                exp.obs.read(
                    value_filter=value_filter or None,
                    column_names=["dataset_id"],
                )
                .concat()
                .to_pandas()
            )
            unique_ids = sorted(df["dataset_id"].unique().tolist())
            logger.info(f"Found {len(unique_ids)} unique dataset(s) "
                        f"matching filters (from {len(df):,} cells)")

            if spec.output.outpath:
                outpath = _resolve_outpath(spec.output.outpath, spec.output.overwrite)
                outpath.write_text("\n".join(unique_ids) + "\n", encoding="utf-8")
                logger.info(f"Wrote dataset list to {outpath}")

            return unique_ids

        # Determine export columns (not needed for dataset_list)
        organ_col = choose_organ_column(obs_keys)
        export_cols = build_obs_export_columns(
            obs_keys=obs_keys,
            export_all_non_ontology=spec.export_all_non_ontology_obs_columns,
            user_cols=spec.obs_columns,
        )

        # Ensure core metadata columns are present
        for required in ("cell_type", "tissue", organ_col):
            if required in obs_keys and required not in export_cols:
                export_cols.append(required)

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

        raise ValueError(f"Unknown output mode: {mode}")


def main():
    """Command-line interface for census_query."""
    parser = argparse.ArgumentParser(
        description="Query CELLxGENE Census using YAML configuration.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Example:
    python -m census_query --config census_query.yaml

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

        if isinstance(result, list):
            print(f"\nFound {len(result)} unique dataset(s):")
            for ds_id in result:
                print(f"  {ds_id}")
        elif isinstance(result, pd.DataFrame):
            print(f"\nQuery returned DataFrame with {len(result):,} rows")
            print(f"Columns: {list(result.columns)}")
            print(f"\nFirst 5 rows:\n{result.head()}")
        else:
            print(f"\nQuery returned: {type(result)}")

    except Exception as e:
        logger.error(f"Query failed: {e}")
        raise SystemExit(1)
