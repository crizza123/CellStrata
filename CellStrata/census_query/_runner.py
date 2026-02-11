"""
_runner.py - Query orchestration and CLI entry point.

This module ties together _schema, _filters, _io, and _visualize to execute
Census queries. It is the only module that imports from all three sub-modules.
"""

from __future__ import annotations

import argparse
import logging
from pathlib import Path
from typing import (
    Any,
    List,
    Union,
)

import pandas as pd
import pyarrow as pa

import cellxgene_census as cxc

from ._schema import QuerySpec, load_query_spec_yaml
from ._filters import (
    build_obs_value_filter,
    build_obs_export_columns,
    choose_organ_column,
)
from ._io import (
    _make_soma_context,
    _resolve_outpath,
    stream_obs_tables,
    write_parquet_stream,
    write_parquet_parts,
)

logger = logging.getLogger(__name__)


def run_query(spec: QuerySpec) -> Union[pd.DataFrame, pa.Table, Any, Path, List[str]]:
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
        - "parquet_dir": Path to output directory of Parquet part files
        - "dataset_list": sorted list of unique dataset ID strings

    Raises:
        ValueError: If output mode is invalid or required options are missing.
    """
    logger.info(f"Starting Census query: version={spec.target.census_version}, "
                f"organism={spec.target.organism}")

    context = _make_soma_context(spec.tiledb_config)

    with cxc.open_soma(census_version=spec.target.census_version, context=context) as census:
        logger.info("Census connection open. Reading experiment schema...")
        exp = census["census_data"][spec.target.organism]
        obs_keys = list(exp.obs.keys())
        logger.info(f"Experiment obs schema: {len(obs_keys)} columns available")

        # Build the value filter
        logger.info("Resolving ontology labels and building value filter...")
        value_filter = build_obs_value_filter(census, spec.obs_filters)
        logger.info(f"Value filter: {value_filter}")

        mode = spec.output.mode
        logger.info(f"Output mode: {mode}")

        # dataset_list mode only needs dataset_id — skip column/organ setup
        if mode == "dataset_list":
            logger.info("Querying Census for dataset IDs (streaming)...")
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

        logger.info(f"Export columns ({len(export_cols)}): {export_cols}")

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

            logger.info("Downloading AnnData from Census (this may take a while)...")
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
            logger.info("Streaming obs metadata from Census...")
            df = (
                exp.obs.read(value_filter=value_filter or None, column_names=export_cols)
                .concat()
                .to_pandas()
            )
            logger.info(f"Retrieved DataFrame: {len(df):,} rows")
            return df

        if mode == "arrow":
            logger.info("Streaming Arrow tables from Census...")
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
            logger.info(f"Streaming obs to Parquet file: {outpath}")

            tables = stream_obs_tables(
                exp, value_filter=value_filter, column_names=export_cols
            )
            rows, batches = write_parquet_stream(
                tables, outpath, compression=spec.output.parquet_compression
            )
            logger.info(f"Wrote Parquet: {rows:,} rows in {batches} batches to {outpath}")
            return outpath

        if mode == "parquet_dir":
            if not spec.output.outpath:
                raise ValueError("output.outpath is required for mode='parquet_dir'")
            outdir = Path(spec.output.outpath)
            logger.info(f"Streaming obs to Parquet directory: {outdir}")

            tables = stream_obs_tables(
                exp, value_filter=value_filter, column_names=export_cols
            )
            rows, parts = write_parquet_parts(
                tables, outdir, compression=spec.output.parquet_compression
            )
            logger.info(
                f"Wrote Parquet directory: {rows:,} rows in {parts} parts to {outdir}"
            )
            return outdir

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
        elif isinstance(result, pa.Table):
            print(f"\nQuery returned Arrow Table with {result.num_rows:,} rows")
        elif isinstance(result, Path):
            print(f"\nQuery output written to: {result}")
        else:
            print(f"\nQuery returned: {type(result)}")

    except Exception as e:
        logger.error(f"Query failed: {e}")
        raise SystemExit(1)
