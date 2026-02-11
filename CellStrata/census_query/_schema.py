"""
_schema.py - Constants, configuration dataclasses, and YAML loading.

This module has no internal dependencies and does not require pandas,
pyarrow, or cellxgene_census, making it lightweight for config-only usage.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from pathlib import Path
from typing import (
    Any,
    Dict,
    List,
    Literal,
    Optional,
    Union,
)

import yaml

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
    mode: Literal["pandas", "arrow", "anndata", "parquet", "dataset_list"] = "pandas"
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
        raw = yaml.safe_load(f) or {}

    target = CensusTarget(**raw.get("target", {}))
    obs_filters = ObsFilters(**raw.get("obs_filters", {}))

    # Filter output keys to those accepted by OutputSpec (ignore removed fields)
    output_raw = raw.get("output", {})
    output_fields = {f.name for f in OutputSpec.__dataclass_fields__.values()}
    output = OutputSpec(**{k: v for k, v in output_raw.items() if k in output_fields})

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
