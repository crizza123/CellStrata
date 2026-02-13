#!/usr/bin/env python3
"""
test_census_query.py - Tests for the CellStrata census package.

Run with: pytest CellStrata/tests/test_census_query.py -v
"""

import sys
from pathlib import Path
from unittest.mock import MagicMock, patch

import pandas as pd
import pytest
import yaml

# Ensure CellStrata package is importable when running tests directly
_repo_root = str(Path(__file__).resolve().parent.parent.parent)
if _repo_root not in sys.path:
    sys.path.insert(0, _repo_root)

# Mock cellxgene_census before importing modules that depend on it,
# so tests can run without the (heavy) Census SDK installed.
sys.modules.setdefault("cellxgene_census", MagicMock())

from CellStrata.census.config import (
    CensusTarget,
    CensusQueryConfig,
    ObsFilters,
    OutputSpec,
    load_config,
)
from CellStrata.census.filters import _q, _in, build_value_filter
from CellStrata.census.query import fetch_obs_df, summarize_by_dataset, validate_obs_cols
from CellStrata.census.outputs import write_summary_csv, write_celllevel_csv


# =========================================================================
# config.py tests
# =========================================================================
class TestCensusTarget:
    def test_defaults(self):
        t = CensusTarget()
        assert t.census_version == "stable"
        assert t.organism == "homo_sapiens"

    def test_custom(self):
        t = CensusTarget(census_version="2024-07-01", organism="mus_musculus")
        assert t.census_version == "2024-07-01"
        assert t.organism == "mus_musculus"

    def test_immutability(self):
        t = CensusTarget()
        with pytest.raises(AttributeError):
            t.census_version = "new"


class TestObsFilters:
    def test_defaults(self):
        f = ObsFilters()
        assert f.is_primary_data is True
        assert f.dataset_id is None
        assert f.sex == []
        assert f.disease == []
        assert f.tissue_general == []
        assert f.cell_type == []
        assert f.assay == []
        assert f.extra_value_filter is None

    def test_custom_lists(self):
        f = ObsFilters(sex=["male"], disease=["normal", "COVID-19"])
        assert f.sex == ["male"]
        assert len(f.disease) == 2


class TestOutputSpec:
    def test_defaults(self):
        o = OutputSpec()
        assert o.summary_csv == "results.csv"
        assert o.celllevel_csv is None
        assert o.write_celllevel is False

    def test_custom(self):
        o = OutputSpec(summary_csv="out.csv", celllevel_csv="cells.csv", write_celllevel=True)
        assert o.summary_csv == "out.csv"
        assert o.write_celllevel is True


class TestCensusQueryConfig:
    def test_defaults(self):
        cfg = CensusQueryConfig()
        assert isinstance(cfg.target, CensusTarget)
        assert isinstance(cfg.obs_filters, ObsFilters)
        assert isinstance(cfg.output, OutputSpec)
        assert "dataset_id" in cfg.obs_cols

    def test_obs_cols_default_has_key_columns(self):
        cfg = CensusQueryConfig()
        for col in ["dataset_id", "cell_type", "tissue_general", "sex"]:
            assert col in cfg.obs_cols


class TestLoadConfig:
    def test_full_config(self, tmp_path):
        cfg_data = {
            "target": {"census_version": "2024-07-01", "organism": "mus_musculus"},
            "obs_filters": {
                "is_primary_data": True,
                "disease": ["normal"],
                "tissue_general": ["brain"],
                "cell_type": ["T cell"],
            },
            "output": {
                "summary_csv": "out.csv",
                "write_celllevel": True,
                "celllevel_csv": "cells.csv",
            },
        }
        p = tmp_path / "cfg.yaml"
        p.write_text(yaml.dump(cfg_data))
        cfg = load_config(p)
        assert cfg.target.organism == "mus_musculus"
        assert cfg.obs_filters.disease == ["normal"]
        assert cfg.output.summary_csv == "out.csv"
        assert cfg.output.write_celllevel is True

    def test_minimal_config(self, tmp_path):
        p = tmp_path / "cfg.yaml"
        p.write_text(yaml.dump({}))
        cfg = load_config(p)
        assert cfg.target.census_version == "stable"
        assert cfg.obs_filters.is_primary_data is True

    def test_empty_file(self, tmp_path):
        """Empty YAML should return all-defaults config."""
        p = tmp_path / "empty.yaml"
        p.write_text("")
        cfg = load_config(p)
        assert cfg.target.census_version == "stable"

    def test_missing_file_raises(self):
        with pytest.raises(FileNotFoundError):
            load_config("/nonexistent/path.yaml")


# =========================================================================
# filters.py tests
# =========================================================================
class TestQ:
    def test_simple_string(self):
        assert _q("hello") == "'hello'"

    def test_escapes_single_quotes(self):
        result = _q("10x 3' v3")
        assert "\\'" in result


class TestIn:
    def test_empty_list_returns_empty(self):
        assert _in("sex", []) == ""

    def test_single_value(self):
        result = _in("sex", ["male"])
        assert result == "sex == 'male'"

    def test_multiple_values(self):
        result = _in("sex", ["male", "female"])
        assert "sex in [" in result
        assert "'male'" in result
        assert "'female'" in result

    def test_none_values_filtered_out(self):
        result = _in("sex", [None, "male", ""])
        assert result == "sex == 'male'"


class TestBuildValueFilter:
    def test_defaults_only(self):
        f = ObsFilters()
        result = build_value_filter(f)
        assert "is_primary_data == True" in result

    def test_with_dataset_id(self):
        f = ObsFilters(dataset_id="abc-123")
        result = build_value_filter(f)
        assert "dataset_id == 'abc-123'" in result

    def test_with_lists(self):
        f = ObsFilters(sex=["male"], disease=["normal"], tissue_general=["brain"])
        result = build_value_filter(f)
        assert "sex == 'male'" in result
        assert "disease == 'normal'" in result
        assert "tissue_general == 'brain'" in result

    def test_with_extra_value_filter(self):
        f = ObsFilters(extra_value_filter="donor_age >= 18")
        result = build_value_filter(f)
        assert "(donor_age >= 18)" in result

    def test_parts_joined_with_and(self):
        f = ObsFilters(sex=["male"], disease=["normal"])
        result = build_value_filter(f)
        assert " and " in result

    def test_multiple_values_use_in(self):
        f = ObsFilters(sex=["male", "female"])
        result = build_value_filter(f)
        assert "sex in [" in result


# =========================================================================
# query.py tests
# =========================================================================
class TestValidateObsCols:
    def test_valid_cols_pass(self):
        mock_census = MagicMock()
        mock_exp = MagicMock()
        mock_exp.obs.keys.return_value = ["dataset_id", "cell_type", "sex"]
        mock_census.__getitem__ = MagicMock(
            side_effect=lambda k: {"census_data": {"homo_sapiens": mock_exp}}[k]
        )
        validate_obs_cols(mock_census, "homo_sapiens", ["dataset_id", "cell_type"])

    def test_missing_cols_raise(self):
        mock_census = MagicMock()
        mock_exp = MagicMock()
        mock_exp.obs.keys.return_value = ["dataset_id"]
        mock_census.__getitem__ = MagicMock(
            side_effect=lambda k: {"census_data": {"homo_sapiens": mock_exp}}[k]
        )
        with pytest.raises(ValueError, match="Missing obs columns"):
            validate_obs_cols(mock_census, "homo_sapiens", ["dataset_id", "nonexistent"])


class TestFetchObsDf:
    @patch("CellStrata.census.query.cxc")
    def test_returns_df_and_filter(self, mock_cxc):
        sample_df = pd.DataFrame({
            "dataset_id": ["d1", "d2"],
            "cell_type": ["T cell", "B cell"],
        })

        mock_exp = MagicMock()
        mock_exp.obs.keys.return_value = [
            "dataset_id", "donor_id", "assay", "sex",
            "tissue", "tissue_general", "cell_type",
            "disease", "development_stage",
        ]
        mock_exp.obs.read.return_value.concat.return_value.to_pandas.return_value = sample_df

        mock_census = MagicMock()
        mock_census.__getitem__ = MagicMock(
            side_effect=lambda k: {"census_data": {"homo_sapiens": mock_exp}}[k]
        )

        mock_cxc.open_soma.return_value.__enter__ = MagicMock(return_value=mock_census)
        mock_cxc.open_soma.return_value.__exit__ = MagicMock(return_value=False)

        cfg = CensusQueryConfig()
        df, vf = fetch_obs_df(cfg)
        assert len(df) == 2
        assert "is_primary_data" in vf


class TestSummarizeByDataset:
    def test_groups_and_counts(self):
        df = pd.DataFrame({"dataset_id": ["a", "a", "b", "b", "b"]})
        summary = summarize_by_dataset(df)
        assert len(summary) == 2
        assert summary.loc[summary["dataset_id"] == "a", "cell_count"].values[0] == 2
        assert summary.loc[summary["dataset_id"] == "b", "cell_count"].values[0] == 3

    def test_empty_df(self):
        df = pd.DataFrame({"dataset_id": pd.Series([], dtype=str)})
        summary = summarize_by_dataset(df)
        assert len(summary) == 0
        assert list(summary.columns) == ["dataset_id", "cell_count"]


# =========================================================================
# outputs.py tests
# =========================================================================
class TestWriteSummaryCsv:
    def test_writes_file(self, tmp_path):
        df = pd.DataFrame({"dataset_id": ["a"], "cell_count": [100]})
        out = tmp_path / "summary.csv"
        write_summary_csv(df, out)
        assert out.exists()
        result = pd.read_csv(out)
        assert result["cell_count"].iloc[0] == 100

    def test_creates_parent_dirs(self, tmp_path):
        out = tmp_path / "sub" / "dir" / "summary.csv"
        df = pd.DataFrame({"dataset_id": ["a"], "cell_count": [10]})
        write_summary_csv(df, out)
        assert out.exists()


class TestWriteCelllevelCsv:
    def test_writes_file(self, tmp_path):
        df = pd.DataFrame({"dataset_id": ["a", "b"], "cell_type": ["T cell", "B cell"]})
        out = tmp_path / "cells.csv"
        write_celllevel_csv(df, out)
        assert out.exists()
        result = pd.read_csv(out)
        assert len(result) == 2
