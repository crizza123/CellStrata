#!/usr/bin/env python3
"""
test_census_query.py - Test suite for the CellStrata census_query module.

#Disclaimer AI was 


This test suite covers:
- Configuration dataclasses
- YAML loading and validation
- Filter building logic
- Helper functions
- Integration tests (requires network access)

Run with: pytest tests/test_census_query.py -v
Run unit tests only: pytest tests/test_census_query.py -v -m "not integration"
"""

import os
import tempfile
from pathlib import Path
from unittest.mock import MagicMock, patch

import pandas as pd
import pyarrow as pa
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
    _resolve_outpath,
    stream_obs_tables,
    write_parquet_stream,
    run_query,
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
    mock = MagicMock()
    
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
    mock["census_info"]["summary_cell_counts"] = mock_scc
    
    # Mock experiment obs keys
    mock_exp = MagicMock()
    mock_exp.obs.keys.return_value = [
        "dataset_id", "donor_id", "assay", "cell_type", "tissue",
        "tissue_general", "disease", "development_stage", "sex",
        "is_primary_data", "suspension_type", "soma_joinid",
        "assay_ontology_term_id", "cell_type_ontology_term_id",
    ]
    mock["census_data"]["homo_sapiens"] = mock_exp
    
    return mock


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
        assert output.parquet_compression == "zstd"
    
    def test_parquet_mode(self):
        """Test parquet output configuration."""
        output = OutputSpec(
            mode="parquet",
            outpath="/tmp/output.parquet",
            overwrite=False,
            parquet_compression="snappy",
        )
        assert output.mode == "parquet"
        assert output.outpath == "/tmp/output.parquet"
        assert output.overwrite is False
        assert output.parquet_compression == "snappy"


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
            output=OutputSpec(mode="arrow"),
        )
        assert spec.target.organism == "mus_musculus"
        assert spec.obs_filters.cell_type_labels == ["mast cell"]
        assert spec.obs_columns == ["dataset_id", "cell_type"]
        assert spec.tiledb_config["vfs.s3.connect_timeout_ms"] == 60000
        assert spec.output.mode == "arrow"


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
            yaml.dump({"output": {"mode": "arrow"}}, f)
        
        spec = load_query_spec_yaml(config_path)
        
        # Check defaults are applied
        assert spec.target.census_version == "stable"
        assert spec.obs_filters.is_primary_data is True
        assert spec.output.mode == "arrow"
    
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
# Unit Tests: I/O Functions
# =============================================================================
class TestResolveOutpath:
    """Tests for _resolve_outpath function."""
    
    def test_creates_parent_directory(self, tmp_path):
        """Test that parent directories are created."""
        outpath = tmp_path / "subdir" / "output.parquet"
        result = _resolve_outpath(str(outpath), overwrite=True)
        
        assert result.parent.exists()
    
    def test_overwrite_existing_file(self, tmp_path):
        """Test overwriting an existing file."""
        outpath = tmp_path / "existing.parquet"
        outpath.touch()
        
        result = _resolve_outpath(str(outpath), overwrite=True)
        assert not outpath.exists()  # File should be deleted
    
    def test_no_overwrite_raises_error(self, tmp_path):
        """Test error when file exists and overwrite=False."""
        outpath = tmp_path / "existing.parquet"
        outpath.touch()
        
        with pytest.raises(FileExistsError):
            _resolve_outpath(str(outpath), overwrite=False)


class TestWriteParquetStream:
    """Tests for write_parquet_stream function."""
    
    def test_write_single_table(self, tmp_path):
        """Test writing a single Arrow table."""
        outpath = tmp_path / "output.parquet"
        table = pa.table({"col1": [1, 2, 3], "col2": ["a", "b", "c"]})
        
        rows, batches = write_parquet_stream([table], outpath)
        
        assert rows == 3
        assert batches == 1
        assert outpath.exists()
    
    def test_write_multiple_tables(self, tmp_path):
        """Test writing multiple Arrow tables."""
        outpath = tmp_path / "output.parquet"
        tables = [
            pa.table({"col1": [1, 2], "col2": ["a", "b"]}),
            pa.table({"col1": [3, 4], "col2": ["c", "d"]}),
        ]
        
        rows, batches = write_parquet_stream(tables, outpath)
        
        assert rows == 4
        assert batches == 2
    
    def test_skip_empty_tables(self, tmp_path):
        """Test that empty tables are skipped."""
        outpath = tmp_path / "output.parquet"
        tables = [
            pa.table({"col1": [], "col2": []}),
            pa.table({"col1": [1], "col2": ["a"]}),
        ]
        
        rows, batches = write_parquet_stream(tables, outpath)
        
        assert rows == 1
        assert batches == 1
    
    def test_compression_option(self, tmp_path):
        """Test different compression options."""
        outpath = tmp_path / "output.parquet"
        table = pa.table({"col1": [1, 2, 3]})
        
        write_parquet_stream([table], outpath, compression="snappy")
        
        assert outpath.exists()


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
    
    @pytest.mark.skip(reason="Requires network access to CELLxGENE Census")
    def test_full_query_arrow(self, integration_config, tmp_path):
        """Test a full query returning Arrow Table."""
        spec = load_query_spec_yaml(integration_config)
        spec = QuerySpec(
            target=spec.target,
            obs_filters=spec.obs_filters,
            output=OutputSpec(mode="arrow"),
        )
        result = run_query(spec)
        
        assert isinstance(result, pa.Table)
        assert result.num_rows > 0


# =============================================================================
# Edge Cases and Error Handling
# =============================================================================
class TestDatasetListMode:
    """Tests for the dataset_list output mode."""

    def test_returns_unique_sorted_ids(self, mock_census):
        """Test that dataset_list mode returns sorted unique dataset IDs."""
        # Mock obs.read to return duplicate dataset_id values
        mock_exp = mock_census["census_data"]["homo_sapiens"]
        obs_df = pd.DataFrame({"dataset_id": ["ds_b", "ds_a", "ds_b", "ds_c", "ds_a"]})
        mock_exp.obs.read.return_value.concat.return_value.to_pandas.return_value = obs_df

        spec = QuerySpec(output=OutputSpec(mode="dataset_list"))

        with patch("census_query.cxc.open_soma") as mock_open:
            mock_open.return_value.__enter__.return_value = mock_census
            result = run_query(spec)

        assert result == ["ds_a", "ds_b", "ds_c"]

    def test_queries_only_dataset_id_column(self, mock_census):
        """Test that only the dataset_id column is requested for efficiency."""
        mock_exp = mock_census["census_data"]["homo_sapiens"]
        obs_df = pd.DataFrame({"dataset_id": ["ds_a"]})
        mock_exp.obs.read.return_value.concat.return_value.to_pandas.return_value = obs_df

        spec = QuerySpec(output=OutputSpec(mode="dataset_list"))

        with patch("census_query.cxc.open_soma") as mock_open:
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

        with patch("census_query.cxc.open_soma") as mock_open:
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

        with patch("census_query.cxc.open_soma") as mock_open:
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
        
        with patch("census_query.cxc.open_soma") as mock_open:
            mock_open.return_value.__enter__.return_value = mock_census
            with pytest.raises(ValueError, match="Unknown output mode"):
                run_query(spec)
    
    def test_parquet_mode_without_outpath(self, mock_census):
        """Test error when parquet mode lacks outpath."""
        spec = QuerySpec(
            output=OutputSpec(mode="parquet", outpath=None),
        )
        
        with patch("census_query.cxc.open_soma") as mock_open:
            mock_open.return_value.__enter__.return_value = mock_census
            with pytest.raises(ValueError, match="outpath is required"):
                run_query(spec)
    
    def test_yaml_with_invalid_syntax(self, tmp_path):
        """Test handling of malformed YAML."""
        config_path = tmp_path / "invalid.yaml"
        with open(config_path, "w") as f:
            f.write("invalid: yaml: content: [")
        
        with pytest.raises(yaml.YAMLError):
            load_query_spec_yaml(config_path)


# =============================================================================
# Test Runner Configuration
# =============================================================================
if __name__ == "__main__":
    # Run tests with pytest
    pytest.main([__file__, "-v", "-m", "not integration"])