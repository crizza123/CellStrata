"""
census_query - YAML-driven CELLxGENE Census query sub-package for CellStrata.

#Disclaimer AI was used to help write this module. Please review carefully.

This package provides functionality to query and download single-cell RNA-seq
data from the CELLxGENE Census database. It supports:
  - YAML configuration for reproducible queries
  - Label-to-ontology-term-id resolution
  - Multiple output formats: pandas DataFrame, Arrow Table, AnnData, Parquet
  - Streaming large queries to avoid OOM errors
  - TileDB context tuning for S3 timeout mitigation

Sub-modules:
    _schema:   Constants, configuration dataclasses, YAML loading
    _filters:  Filter building and ontology term resolution
    _io:       I/O utilities, streaming, file management
    _runner:   Query orchestration and CLI entry point

Usage:
    from census_query import load_query_spec_yaml, run_query
    spec = load_query_spec_yaml("config/census_query.yaml")
    result = run_query(spec)

Author: Crizza
Course: BIOINF 576
"""

import logging

# Configure logging once at package level
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
)

# --- Schema: constants, dataclasses, YAML loading ---
from ._schema import (
    OBS_TERM_ID_COLS,
    OBS_LABEL_COLS,
    CensusTarget,
    ObsFilters,
    OutputSpec,
    QuerySpec,
    load_query_spec_yaml,
)

# --- Filters: filter building and ontology resolution ---
from ._filters import (
    _format_in_list,
    _get_summary_cell_counts,
    resolve_ontology_term_ids,
    choose_organ_column,
    build_obs_value_filter,
    _default_obs_columns,
    build_obs_export_columns,
)

# --- I/O: streaming, file management, TileDB context ---
from ._io import (
    _resolve_outpath,
    stream_obs_tables,
    write_parquet_stream,
)

# --- Runner: query orchestration and CLI ---
from ._runner import run_query, main

__all__ = [
    # Dataclasses
    "QuerySpec",
    "CensusTarget",
    "ObsFilters",
    "OutputSpec",
    # Schema constants
    "OBS_TERM_ID_COLS",
    "OBS_LABEL_COLS",
    # Config loading
    "load_query_spec_yaml",
    # Filters
    "resolve_ontology_term_ids",
    "choose_organ_column",
    "build_obs_value_filter",
    "build_obs_export_columns",
    # I/O
    "stream_obs_tables",
    "write_parquet_stream",
    # Runner
    "run_query",
    "main",
]
