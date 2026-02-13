from __future__ import annotations
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional
import yaml


@dataclass(frozen=True)
class CensusTarget:
    census_version: str = "stable"
    organism: str = "homo_sapiens"


@dataclass(frozen=True)
class ObsFilters:
    is_primary_data: bool = True
    dataset_id: Optional[str] = None

    sex: List[str] = field(default_factory=list)
    disease: List[str] = field(default_factory=list)          # prefer label like "normal"
    tissue_general: List[str] = field(default_factory=list)
    cell_type: List[str] = field(default_factory=list)
    assay: List[str] = field(default_factory=list)

    extra_value_filter: Optional[str] = None                  # raw appended clause


@dataclass(frozen=True)
class OutputSpec:
    summary_csv: str = "results.csv"
    celllevel_csv: Optional[str] = None
    write_celllevel: bool = False


@dataclass(frozen=True)
class CensusQueryConfig:
    target: CensusTarget = CensusTarget()
    obs_filters: ObsFilters = ObsFilters()
    output: OutputSpec = OutputSpec()
    obs_cols: List[str] = field(default_factory=lambda: [
        "dataset_id", "donor_id", "assay", "sex",
        "tissue", "tissue_general", "cell_type",
        "disease", "development_stage",
    ])


def load_config(path: str | Path) -> CensusQueryConfig:
    raw = yaml.safe_load(Path(path).read_text())
    cfg = raw if raw is not None else {}

    target = CensusTarget(**cfg.get("target", {}))
    obs_filters = ObsFilters(**cfg.get("obs_filters", {}))
    output = OutputSpec(**cfg.get("output", {}))
    obs_cols = cfg.get("obs_cols") or CensusQueryConfig().obs_cols

    return CensusQueryConfig(
        target=target,
        obs_filters=obs_filters,
        output=output,
        obs_cols=obs_cols,
    )
