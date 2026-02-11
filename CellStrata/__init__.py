"""
CellStrata - Single-cell RNA-seq data acquisition and analysis toolkit.

Disclaimer this module was assisted by AI. Please review the code carefully.

This package provides tools for:
- Querying and downloading data from CELLxGENE Census
- Visualizing cell metadata distributions
- Quality control and normalization
- Cell type annotation with focus on mast cell heterogeneity

Modules:
    census_query: YAML-driven CELLxGENE Census query functionality
"""

__version__ = "0.1.0"
__author__ = "Crizza"

from .census_query import (
    QuerySpec,
    CensusTarget,
    ObsFilters,
    OutputSpec,
    load_query_spec_yaml,
    run_query,
    plot_cell_type_counts,
    plot_tissue_counts,
    plot_sex_distribution,
    plot_assay_counts,
    plot_disease_counts,
    plot_dataset_contribution,
    plot_development_stage_counts,
    plot_donors_per_dataset,
    plot_cell_type_by_tissue,
    plot_metadata_summary,
)

__all__ = [
    "QuerySpec",
    "CensusTarget",
    "ObsFilters",
    "OutputSpec",
    "load_query_spec_yaml",
    "run_query",
    "plot_cell_type_counts",
    "plot_tissue_counts",
    "plot_sex_distribution",
    "plot_assay_counts",
    "plot_disease_counts",
    "plot_dataset_contribution",
    "plot_development_stage_counts",
    "plot_donors_per_dataset",
    "plot_cell_type_by_tissue",
    "plot_metadata_summary",
]
