#!/usr/bin/env python3
"""
test_census_query.py - Test suite for the CellStrata census_query module.

#Disclaimer AI was used to help write this module. Please review carefully.

This test suite covers:
- Configuration dataclasses
- YAML loading and validation
- Filter building logic
- Helper functions
- Visualization functions
- Integration tests (requires network access)

Run with: pytest tests/test_census_query.py -v
Run unit tests only: pytest tests/test_census_query.py -v -m "not integration"
"""

import os
import tempfile
from pathlib import Path
from unittest.mock import MagicMock, patch

import matplotlib
matplotlib.use("Agg")  # non-interactive backend for CI
import matplotlib.pyplot as plt

import pandas as pd
import pytest
import yaml

# Import module under test
import sys
sys.path.insert(0, str(Path(__file__).parent.parent))

from census_query import (
    CensusTarget,
    ObsFilters,
    OutputSpec,
    QuerySpec,
    load_query_spec_yaml,
    _format_in_list,
    _default_obs_columns,
    build_obs_export_columns,
    choose_organ_column,
    build_obs_value_filter,
    resolve_ontology_term_ids,
    run_query,
    # Visualization functions
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


# =============================================================================
# Fixtures
# =============================================================================
@pytest.fixture
def sample_yaml_config():
    """Create a sample YAML configuration for testing."""
    return {
        "target": {
            "census_version": "stable",
            "organism": "homo_sapiens",
            "measurement": "RNA",
        },
        "obs_filters": {
            "is_primary_data": True,
            "suspension_type": "cell",
            "disease_ontology_term_ids": ["PATO:0000461"],
            "assay_labels": ["10x 3' v3"],
            "tissue_general_labels": ["lung"],
            "cell_type_labels": ["B cell"],
        },
        "output": {
            "mode": "pandas",
        },
    }


@pytest.fixture
def yaml_config_file(sample_yaml_config, tmp_path):
    """Create a temporary YAML config file."""
    config_path = tmp_path / "test_config.yaml"
    with open(config_path, "w") as f:
        yaml.dump(sample_yaml_config, f)
    return config_path


@pytest.fixture
def mock_census():
    """Create a mock Census object for testing."""
    # Mock summary_cell_counts table
    scc_data = pd.DataFrame({
        "category": ["disease", "disease", "assay", "assay", "tissue_general",
                     "cell_type", "sex", "sex"],
        "label": ["normal", "COVID-19", "10x 3' v3", "10x 3' v2", "lung",
                  "B cell", "male", "female"],
        "ontology_term_id": ["PATO:0000461", "MONDO:0100096", "EFO:0009922",
                            "EFO:0009899", "UBERON:0002048", "CL:0000236",
                            "PATO:0000384", "PATO:0000383"],
    })

    mock_scc = MagicMock()
    mock_scc.read.return_value.concat.return_value.to_pandas.return_value = scc_data

    # Mock experiment obs keys
    mock_exp = MagicMock()
    mock_exp.obs.keys.return_value = [
        "dataset_id", "donor_id", "assay", "cell_type", "tissue",
        "tissue_general", "disease", "development_stage", "sex",
        "is_primary_data", "suspension_type", "soma_joinid",
        "assay_ontology_term_id", "cell_type_ontology_term_id",
    ]

    # Build census mock with explicit __getitem__ routing
    census_info = MagicMock()
    census_info.__getitem__ = MagicMock(
        side_effect=lambda key: {"summary_cell_counts": mock_scc}[key]
    )

    census_data = MagicMock()
    census_data.__getitem__ = MagicMock(
        side_effect=lambda key: {"homo_sapiens": mock_exp}[key]
    )

    mock = MagicMock()
    mock.__getitem__ = MagicMock(
        side_effect=lambda key: {
            "census_info": census_info,
            "census_data": census_data,
        }[key]
    )

    return mock


@pytest.fixture
def sample_metadata_df():
    """Create a sample metadata DataFrame for visualization tests."""
    return pd.DataFrame({
        "dataset_id": (["ds_a"] * 80 + ["ds_b"] * 60 + ["ds_c"] * 40 + ["ds_d"] * 20),
        "donor_id": (["donor_1"] * 50 + ["donor_2"] * 50 +
                     ["donor_3"] * 50 + ["donor_4"] * 50),
        "cell_type": (["B cell"] * 60 + ["T cell"] * 50 + ["mast cell"] * 40 +
                      ["macrophage"] * 30 + ["NK cell"] * 20),
        "tissue_general": (["blood"] * 70 + ["lung"] * 60 +
                           ["skin"] * 40 + ["heart"] * 30),
        "tissue": (["peripheral blood"] * 70 + ["lung parenchyma"] * 60 +
                   ["skin of body"] * 40 + ["heart left ventricle"] * 30),
        "disease": ["normal"] * 150 + ["COVID-19"] * 50,
        "sex": ["male"] * 100 + ["female"] * 100,
        "assay": ["10x 3' v3"] * 120 + ["10x 3' v2"] * 80,
        "development_stage": ["adult"] * 160 + ["child"] * 40,
        "self_reported_ethnicity": ["European"] * 100 + ["Asian"] * 100,
        "is_primary_data": [True] * 200,
        "suspension_type": ["cell"] * 200,
    })


# =============================================================================
# Unit Tests: Dataclasses
# =============================================================================
class TestCensusTarget:
    """Tests for CensusTarget dataclass."""

    def test_default_values(self):
        """Test default initialization."""
        target = CensusTarget()
        assert target.census_version == "stable"
        assert target.organism == "homo_sapiens"
        assert target.measurement == "RNA"

    def test_custom_values(self):
        """Test custom initialization."""
        target = CensusTarget(
            census_version="2025-11-08",
            organism="mus_musculus",
            measurement="ATAC",
        )
        assert target.census_version == "2025-11-08"
        assert target.organism == "mus_musculus"
        assert target.measurement == "ATAC"

    def test_immutability(self):
        """Test that dataclass is frozen/immutable."""
        target = CensusTarget()
        with pytest.raises(AttributeError):
            target.census_version = "new_version"


class TestObsFilters:
    """Tests for ObsFilters dataclass."""

    def test_default_values(self):
        """Test default initialization."""
        filters = ObsFilters()
        assert filters.is_primary_data is True
        assert filters.suspension_type == "cell"
        assert filters.disease_ontology_term_ids == []
        assert filters.cell_type_labels == []

    def test_custom_filters(self):
        """Test custom filter initialization."""
        filters = ObsFilters(
            is_primary_data=False,
            disease_ontology_term_ids=["PATO:0000461"],
            cell_type_labels=["B cell", "T cell"],
        )
        assert filters.is_primary_data is False
        assert filters.disease_ontology_term_ids == ["PATO:0000461"]
        assert filters.cell_type_labels == ["B cell", "T cell"]

    def test_extra_value_filter(self):
        """Test extra_value_filter field."""
        filters = ObsFilters(extra_value_filter="donor_age >= 18")
        assert filters.extra_value_filter == "donor_age >= 18"


class TestOutputSpec:
    """Tests for OutputSpec dataclass."""

    def test_default_values(self):
        """Test default initialization."""
        output = OutputSpec()
        assert output.mode == "pandas"
        assert output.outpath is None
        assert output.overwrite is True

    def test_dataset_list_mode(self):
        """Test dataset_list output configuration."""
        output = OutputSpec(
            mode="dataset_list",
            outpath="/tmp/datasets.txt",
            overwrite=False,
        )
        assert output.mode == "dataset_list"
        assert output.outpath == "/tmp/datasets.txt"
        assert output.overwrite is False


class TestQuerySpec:
    """Tests for QuerySpec dataclass."""

    def test_default_values(self):
        """Test default initialization."""
        spec = QuerySpec()
        assert isinstance(spec.target, CensusTarget)
        assert isinstance(spec.obs_filters, ObsFilters)
        assert isinstance(spec.output, OutputSpec)
        assert spec.obs_columns == []
        assert spec.tiledb_config == {}

    def test_full_specification(self):
        """Test full query specification."""
        spec = QuerySpec(
            target=CensusTarget(organism="mus_musculus"),
            obs_filters=ObsFilters(cell_type_labels=["mast cell"]),
            obs_columns=["dataset_id", "cell_type"],
            tiledb_config={"vfs.s3.connect_timeout_ms": 60000},
            output=OutputSpec(mode="pandas"),
        )
        assert spec.target.organism == "mus_musculus"
        assert spec.obs_filters.cell_type_labels == ["mast cell"]
        assert spec.obs_columns == ["dataset_id", "cell_type"]
        assert spec.tiledb_config["vfs.s3.connect_timeout_ms"] == 60000
        assert spec.output.mode == "pandas"


# =============================================================================
# Unit Tests: YAML Loading
# =============================================================================
class TestLoadQuerySpecYaml:
    """Tests for YAML configuration loading."""

    def test_load_valid_config(self, yaml_config_file):
        """Test loading a valid YAML configuration."""
        spec = load_query_spec_yaml(yaml_config_file)

        assert spec.target.census_version == "stable"
        assert spec.target.organism == "homo_sapiens"
        assert spec.obs_filters.is_primary_data is True
        assert spec.obs_filters.disease_ontology_term_ids == ["PATO:0000461"]
        assert spec.output.mode == "pandas"

    def test_load_missing_file(self):
        """Test error handling for missing config file."""
        with pytest.raises(FileNotFoundError):
            load_query_spec_yaml("/nonexistent/path/config.yaml")

    def test_load_minimal_config(self, tmp_path):
        """Test loading a minimal configuration with defaults."""
        config_path = tmp_path / "minimal.yaml"
        with open(config_path, "w") as f:
            yaml.dump({"output": {"mode": "pandas"}}, f)

        spec = load_query_spec_yaml(config_path)

        # Check defaults are applied
        assert spec.target.census_version == "stable"
        assert spec.obs_filters.is_primary_data is True
        assert spec.output.mode == "pandas"

    def test_load_empty_config(self, tmp_path):
        """Test loading an empty configuration file."""
        config_path = tmp_path / "empty.yaml"
        config_path.touch()

        # Empty YAML returns None, which should use all defaults
        spec = load_query_spec_yaml(config_path)
        assert spec.target.census_version == "stable"

    def test_load_config_with_tiledb(self, tmp_path):
        """Test loading configuration with TileDB settings."""
        config = {
            "tiledb_config": {
                "vfs.s3.connect_timeout_ms": 60000,
                "vfs.s3.request_timeout_ms": 600000,
            }
        }
        config_path = tmp_path / "tiledb.yaml"
        with open(config_path, "w") as f:
            yaml.dump(config, f)

        spec = load_query_spec_yaml(config_path)
        assert spec.tiledb_config["vfs.s3.connect_timeout_ms"] == 60000


# =============================================================================
# Unit Tests: Helper Functions
# =============================================================================
class TestFormatInList:
    """Tests for _format_in_list helper."""

    def test_simple_list(self):
        """Test formatting a simple list."""
        result = _format_in_list(["a", "b", "c"])
        assert result == "['a', 'b', 'c']"

    def test_single_item(self):
        """Test formatting a single-item list."""
        result = _format_in_list(["only_one"])
        assert result == "['only_one']"

    def test_empty_list(self):
        """Test formatting an empty list."""
        result = _format_in_list([])
        assert result == "[]"

    def test_escape_quotes(self):
        """Test escaping single quotes in values."""
        result = _format_in_list(["10x 3' v3"])
        assert "\\'" in result

    def test_escape_backslash(self):
        """Test escaping backslashes in values."""
        result = _format_in_list(["path\\to\\file"])
        assert "\\\\" in result


class TestDefaultObsColumns:
    """Tests for _default_obs_columns helper."""

    def test_returns_list(self):
        """Test that function returns a list."""
        cols = _default_obs_columns()
        assert isinstance(cols, list)

    def test_contains_expected_columns(self):
        """Test that default columns include expected fields."""
        cols = _default_obs_columns()
        expected = ["dataset_id", "donor_id", "cell_type", "tissue", "disease"]
        for col in expected:
            assert col in cols

    def test_excludes_ontology_ids(self):
        """Test that ontology term ID columns are excluded by default."""
        cols = _default_obs_columns()
        for col in cols:
            assert "ontology_term_id" not in col


class TestChooseOrganColumn:
    """Tests for choose_organ_column helper."""

    def test_prefers_tissue_general(self):
        """Test that tissue_general is preferred over tissue."""
        keys = ["tissue", "tissue_general", "cell_type"]
        result = choose_organ_column(keys)
        assert result == "tissue_general"

    def test_falls_back_to_tissue(self):
        """Test fallback to tissue when tissue_general is absent."""
        keys = ["tissue", "cell_type", "disease"]
        result = choose_organ_column(keys)
        assert result == "tissue"

    def test_raises_when_neither_present(self):
        """Test error when neither tissue column is present."""
        keys = ["cell_type", "disease", "assay"]
        with pytest.raises(ValueError, match="Neither tissue_general nor tissue"):
            choose_organ_column(keys)


class TestBuildObsExportColumns:
    """Tests for build_obs_export_columns helper."""

    def test_user_specified_columns(self):
        """Test using user-specified columns."""
        obs_keys = ["dataset_id", "cell_type", "tissue", "extra_col"]
        result = build_obs_export_columns(
            obs_keys=obs_keys,
            export_all_non_ontology=False,
            user_cols=["dataset_id", "cell_type"],
        )
        assert result == ["dataset_id", "cell_type"]

    def test_export_all_non_ontology(self):
        """Test exporting all non-ontology columns."""
        obs_keys = [
            "dataset_id", "cell_type", "cell_type_ontology_term_id",
            "soma_joinid", "tissue",
        ]
        result = build_obs_export_columns(
            obs_keys=obs_keys,
            export_all_non_ontology=True,
            user_cols=[],
        )
        assert "dataset_id" in result
        assert "cell_type" in result
        assert "tissue" in result
        assert "cell_type_ontology_term_id" not in result
        assert "soma_joinid" not in result

    def test_missing_column_raises_error(self):
        """Test error when requested column is not in schema."""
        obs_keys = ["dataset_id", "cell_type"]
        with pytest.raises(ValueError, match="not present in schema"):
            build_obs_export_columns(
                obs_keys=obs_keys,
                export_all_non_ontology=False,
                user_cols=["nonexistent_column"],
            )


# =============================================================================
# Unit Tests: Filter Building
# =============================================================================
class TestResolveOntologyTermIds:
    """Tests for resolve_ontology_term_ids function."""

    def test_resolve_single_label(self, mock_census):
        """Test resolving a single label."""
        result = resolve_ontology_term_ids(mock_census, "disease", ["normal"])
        assert result == ["PATO:0000461"]

    def test_resolve_multiple_labels(self, mock_census):
        """Test resolving multiple labels."""
        result = resolve_ontology_term_ids(mock_census, "sex", ["male", "female"])
        assert result == ["PATO:0000384", "PATO:0000383"]

    def test_empty_labels(self, mock_census):
        """Test with empty label list."""
        result = resolve_ontology_term_ids(mock_census, "disease", [])
        assert result == []

    def test_missing_label_raises_error(self, mock_census):
        """Test error when label cannot be resolved."""
        with pytest.raises(ValueError, match="Could not resolve"):
            resolve_ontology_term_ids(mock_census, "disease", ["nonexistent_disease"])


class TestBuildObsValueFilter:
    """Tests for build_obs_value_filter function."""

    def test_basic_filters(self, mock_census):
        """Test building filter with basic options."""
        filters = ObsFilters(
            is_primary_data=True,
            suspension_type="cell",
        )
        result = build_obs_value_filter(mock_census, filters)

        assert "is_primary_data == True" in result
        assert "suspension_type == 'cell'" in result

    def test_disease_ontology_filter(self, mock_census):
        """Test disease filter with ontology term IDs."""
        filters = ObsFilters(
            disease_ontology_term_ids=["PATO:0000461"],
        )
        result = build_obs_value_filter(mock_census, filters)

        assert "disease_ontology_term_id in ['PATO:0000461']" in result

    def test_disease_label_resolution(self, mock_census):
        """Test disease filter with label resolution."""
        filters = ObsFilters(
            disease_labels=["normal"],
        )
        result = build_obs_value_filter(mock_census, filters)

        assert "disease_ontology_term_id in ['PATO:0000461']" in result

    def test_multiple_filters_combined(self, mock_census):
        """Test combining multiple filters with 'and'."""
        filters = ObsFilters(
            is_primary_data=True,
            disease_ontology_term_ids=["PATO:0000461"],
            cell_type_labels=["B cell"],
        )
        result = build_obs_value_filter(mock_census, filters)

        assert " and " in result
        assert "is_primary_data == True" in result
        assert "disease_ontology_term_id" in result
        assert "cell_type_ontology_term_id" in result

    def test_extra_value_filter(self, mock_census):
        """Test adding extra custom filter."""
        filters = ObsFilters(
            extra_value_filter="donor_age >= 18",
        )
        result = build_obs_value_filter(mock_census, filters)

        assert "(donor_age >= 18)" in result

    def test_empty_filters(self, mock_census):
        """Test with no filters specified."""
        filters = ObsFilters(
            is_primary_data=None,
            suspension_type=None,
        )
        result = build_obs_value_filter(mock_census, filters)

        assert result == ""


# =============================================================================
# Unit Tests: Runner
# =============================================================================
class TestDatasetListMode:
    """Tests for the dataset_list output mode."""

    def test_returns_unique_sorted_ids(self, mock_census):
        """Test that dataset_list mode returns sorted unique dataset IDs."""
        mock_exp = mock_census["census_data"]["homo_sapiens"]
        obs_df = pd.DataFrame({"dataset_id": ["ds_b", "ds_a", "ds_b", "ds_c", "ds_a"]})
        mock_exp.obs.read.return_value.concat.return_value.to_pandas.return_value = obs_df

        spec = QuerySpec(output=OutputSpec(mode="dataset_list"))

        with patch("census_query._runner.cxc.open_soma") as mock_open:
            mock_open.return_value.__enter__.return_value = mock_census
            result = run_query(spec)

        assert result == ["ds_a", "ds_b", "ds_c"]

    def test_queries_only_dataset_id_column(self, mock_census):
        """Test that only the dataset_id column is requested for efficiency."""
        mock_exp = mock_census["census_data"]["homo_sapiens"]
        obs_df = pd.DataFrame({"dataset_id": ["ds_a"]})
        mock_exp.obs.read.return_value.concat.return_value.to_pandas.return_value = obs_df

        spec = QuerySpec(output=OutputSpec(mode="dataset_list"))

        with patch("census_query._runner.cxc.open_soma") as mock_open:
            mock_open.return_value.__enter__.return_value = mock_census
            run_query(spec)

        mock_exp.obs.read.assert_called_once()
        call_kwargs = mock_exp.obs.read.call_args
        assert call_kwargs.kwargs.get("column_names") == ["dataset_id"]

    def test_writes_to_file_when_outpath_set(self, mock_census, tmp_path):
        """Test that dataset list is written to a file when outpath is provided."""
        mock_exp = mock_census["census_data"]["homo_sapiens"]
        obs_df = pd.DataFrame({"dataset_id": ["ds_b", "ds_a", "ds_b"]})
        mock_exp.obs.read.return_value.concat.return_value.to_pandas.return_value = obs_df

        outpath = tmp_path / "datasets.txt"
        spec = QuerySpec(
            output=OutputSpec(mode="dataset_list", outpath=str(outpath)),
        )

        with patch("census_query._runner.cxc.open_soma") as mock_open:
            mock_open.return_value.__enter__.return_value = mock_census
            result = run_query(spec)

        assert result == ["ds_a", "ds_b"]
        assert outpath.exists()
        lines = outpath.read_text().strip().split("\n")
        assert lines == ["ds_a", "ds_b"]

    def test_returns_empty_list_when_no_matches(self, mock_census):
        """Test that an empty list is returned when no cells match."""
        mock_exp = mock_census["census_data"]["homo_sapiens"]
        obs_df = pd.DataFrame({"dataset_id": pd.Series([], dtype=str)})
        mock_exp.obs.read.return_value.concat.return_value.to_pandas.return_value = obs_df

        spec = QuerySpec(output=OutputSpec(mode="dataset_list"))

        with patch("census_query._runner.cxc.open_soma") as mock_open:
            mock_open.return_value.__enter__.return_value = mock_census
            result = run_query(spec)

        assert result == []


class TestEdgeCases:
    """Tests for edge cases and error handling."""

    def test_invalid_output_mode(self, mock_census):
        """Test error for invalid output mode."""
        spec = QuerySpec(
            output=OutputSpec(mode="invalid_mode"),  # type: ignore
        )

        with patch("census_query._runner.cxc.open_soma") as mock_open:
            mock_open.return_value.__enter__.return_value = mock_census
            with pytest.raises(ValueError):
                run_query(spec)

    def test_yaml_with_invalid_syntax(self, tmp_path):
        """Test handling of malformed YAML."""
        config_path = tmp_path / "invalid.yaml"
        with open(config_path, "w") as f:
            f.write("invalid: yaml: content: [")

        with pytest.raises(yaml.YAMLError):
            load_query_spec_yaml(config_path)


# =============================================================================
# Unit Tests: Visualization
# =============================================================================
class TestPlotCellTypeCounts:
    """Tests for plot_cell_type_counts."""

    def test_returns_axes(self, sample_metadata_df):
        ax = plot_cell_type_counts(sample_metadata_df)
        assert isinstance(ax, plt.Axes)
        plt.close("all")

    def test_respects_top_n(self, sample_metadata_df):
        ax = plot_cell_type_counts(sample_metadata_df, top_n=3)
        # 3 bars in the plot
        assert len(ax.patches) == 3
        plt.close("all")

    def test_accepts_existing_axes(self, sample_metadata_df):
        fig, ax = plt.subplots()
        returned = plot_cell_type_counts(sample_metadata_df, ax=ax)
        assert returned is ax
        plt.close("all")

    def test_missing_column_raises(self):
        df = pd.DataFrame({"tissue": ["lung"]})
        with pytest.raises(KeyError, match="cell_type"):
            plot_cell_type_counts(df)
        plt.close("all")


class TestPlotTissueCounts:
    """Tests for plot_tissue_counts."""

    def test_returns_axes(self, sample_metadata_df):
        ax = plot_tissue_counts(sample_metadata_df)
        assert isinstance(ax, plt.Axes)
        plt.close("all")

    def test_custom_column(self, sample_metadata_df):
        ax = plot_tissue_counts(sample_metadata_df, column="tissue")
        assert isinstance(ax, plt.Axes)
        plt.close("all")

    def test_missing_column_raises(self):
        df = pd.DataFrame({"cell_type": ["B cell"]})
        with pytest.raises(KeyError, match="tissue_general"):
            plot_tissue_counts(df)
        plt.close("all")


class TestPlotSexDistribution:
    """Tests for plot_sex_distribution."""

    def test_returns_axes(self, sample_metadata_df):
        ax = plot_sex_distribution(sample_metadata_df)
        assert isinstance(ax, plt.Axes)
        plt.close("all")

    def test_bar_count_matches_unique_values(self, sample_metadata_df):
        ax = plot_sex_distribution(sample_metadata_df)
        n_unique = sample_metadata_df["sex"].nunique()
        assert len(ax.patches) == n_unique
        plt.close("all")


class TestPlotAssayCounts:
    """Tests for plot_assay_counts."""

    def test_returns_axes(self, sample_metadata_df):
        ax = plot_assay_counts(sample_metadata_df)
        assert isinstance(ax, plt.Axes)
        plt.close("all")


class TestPlotDiseaseCounts:
    """Tests for plot_disease_counts."""

    def test_returns_axes(self, sample_metadata_df):
        ax = plot_disease_counts(sample_metadata_df)
        assert isinstance(ax, plt.Axes)
        plt.close("all")


class TestPlotDatasetContribution:
    """Tests for plot_dataset_contribution."""

    def test_returns_axes(self, sample_metadata_df):
        ax = plot_dataset_contribution(sample_metadata_df)
        assert isinstance(ax, plt.Axes)
        plt.close("all")

    def test_respects_top_n(self, sample_metadata_df):
        ax = plot_dataset_contribution(sample_metadata_df, top_n=2)
        assert len(ax.patches) == 2
        plt.close("all")


class TestPlotDevelopmentStageCounts:
    """Tests for plot_development_stage_counts."""

    def test_returns_axes(self, sample_metadata_df):
        ax = plot_development_stage_counts(sample_metadata_df)
        assert isinstance(ax, plt.Axes)
        plt.close("all")


class TestPlotDonorsPerDataset:
    """Tests for plot_donors_per_dataset."""

    def test_returns_axes(self, sample_metadata_df):
        ax = plot_donors_per_dataset(sample_metadata_df)
        assert isinstance(ax, plt.Axes)
        plt.close("all")


class TestPlotCellTypeByTissue:
    """Tests for plot_cell_type_by_tissue."""

    def test_returns_axes(self, sample_metadata_df):
        ax = plot_cell_type_by_tissue(sample_metadata_df)
        assert isinstance(ax, plt.Axes)
        plt.close("all")

    def test_custom_tissue_column(self, sample_metadata_df):
        ax = plot_cell_type_by_tissue(sample_metadata_df, tissue_column="tissue")
        assert isinstance(ax, plt.Axes)
        plt.close("all")


class TestPlotMetadataSummary:
    """Tests for plot_metadata_summary dashboard."""

    def test_returns_figure(self, sample_metadata_df):
        fig = plot_metadata_summary(sample_metadata_df)
        assert isinstance(fig, plt.Figure)
        plt.close("all")

    def test_save_to_file(self, sample_metadata_df, tmp_path):
        outpath = tmp_path / "summary.png"
        fig = plot_metadata_summary(sample_metadata_df, save_path=outpath)
        assert outpath.exists()
        assert outpath.stat().st_size > 0
        plt.close("all")

    def test_handles_missing_columns_gracefully(self):
        """Dashboard should hide panels for missing columns."""
        df = pd.DataFrame({
            "cell_type": ["B cell", "T cell", "B cell"],
            "tissue_general": ["lung", "blood", "lung"],
        })
        fig = plot_metadata_summary(df)
        assert isinstance(fig, plt.Figure)
        plt.close("all")


# =============================================================================
# Integration Tests (require network access)
# =============================================================================
@pytest.mark.integration
class TestIntegration:
    """Integration tests that require network access to CELLxGENE Census."""

    @pytest.fixture
    def integration_config(self, tmp_path):
        """Create a minimal integration test config."""
        config = {
            "target": {
                "census_version": "stable",
                "organism": "homo_sapiens",
            },
            "obs_filters": {
                "is_primary_data": True,
                "disease_ontology_term_ids": ["PATO:0000461"],
                "cell_type_labels": ["B cell"],
                "tissue_general_labels": ["blood"],
            },
            "output": {
                "mode": "pandas",
            },
        }
        config_path = tmp_path / "integration_config.yaml"
        with open(config_path, "w") as f:
            yaml.dump(config, f)
        return config_path

    @pytest.mark.skip(reason="Requires network access to CELLxGENE Census")
    def test_full_query_pandas(self, integration_config):
        """Test a full query returning pandas DataFrame."""
        spec = load_query_spec_yaml(integration_config)
        result = run_query(spec)

        assert isinstance(result, pd.DataFrame)
        assert len(result) > 0
        assert "cell_type" in result.columns


# =============================================================================
# Test Runner Configuration
# =============================================================================
if __name__ == "__main__":
    # Run tests with pytest
    pytest.main([__file__, "-v", "-m", "not integration"])
