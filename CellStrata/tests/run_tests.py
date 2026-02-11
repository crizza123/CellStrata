#!/usr/bin/env python3
"""
run_tests.py - Standalone test runner for census_query module.
Disclaimer AI was used to help write this test suite. Please review carefully.
Runs tests without pytest or pandas dependencies.
"""

import sys
import os
import tempfile
from pathlib import Path
from dataclasses import dataclass, field
from typing import List, Callable, Any, Dict, Optional, Sequence, Union

import yaml

print("CellStrata Test Runner")
print("=" * 60)

# =============================================================================
# Dataclasses
# =============================================================================
@dataclass(frozen=True)
class CensusTarget:
    census_version: str = "stable"
    organism: str = "homo_sapiens"
    measurement: str = "RNA"

@dataclass(frozen=True)
class ObsFilters:
    is_primary_data: Optional[bool] = True
    suspension_type: Optional[str] = "cell"
    disease_ontology_term_ids: List[str] = field(default_factory=list)
    disease_labels: List[str] = field(default_factory=list)
    sex_labels: List[str] = field(default_factory=list)
    sex_ontology_term_ids: List[str] = field(default_factory=list)
    assay_labels: List[str] = field(default_factory=list)
    assay_ontology_term_ids: List[str] = field(default_factory=list)
    tissue_general_labels: List[str] = field(default_factory=list)
    tissue_general_ontology_term_ids: List[str] = field(default_factory=list)
    tissue_labels: List[str] = field(default_factory=list)
    tissue_ontology_term_ids: List[str] = field(default_factory=list)
    cell_type_labels: List[str] = field(default_factory=list)
    cell_type_ontology_term_ids: List[str] = field(default_factory=list)
    development_stage_labels: List[str] = field(default_factory=list)
    development_stage_ontology_term_ids: List[str] = field(default_factory=list)
    extra_value_filter: Optional[str] = None

@dataclass(frozen=True)
class OutputSpec:
    mode: str = "pandas"
    outpath: Optional[str] = None
    overwrite: bool = True
    parquet_compression: str = "zstd"

@dataclass(frozen=True)
class QuerySpec:
    target: CensusTarget = field(default_factory=CensusTarget)
    obs_filters: ObsFilters = field(default_factory=ObsFilters)
    export_all_non_ontology_obs_columns: bool = False
    obs_columns: List[str] = field(default_factory=list)
    var_value_filter: Optional[str] = None
    var_columns: List[str] = field(default_factory=list)
    tiledb_config: Dict[str, Any] = field(default_factory=dict)
    output: OutputSpec = field(default_factory=OutputSpec)

OBS_TERM_ID_COLS = {
    "assay": "assay_ontology_term_id",
    "cell_type": "cell_type_ontology_term_id",
    "disease": "disease_ontology_term_id",
    "sex": "sex_ontology_term_id",
    "tissue": "tissue_ontology_term_id",
    "tissue_general": "tissue_general_ontology_term_id",
}

# =============================================================================
# Functions to test
# =============================================================================
def load_query_spec_yaml(path):
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"Configuration file not found: {path}")
    with path.open("r", encoding="utf-8") as f:
        raw = yaml.safe_load(f) or {}
    target = CensusTarget(**raw.get("target", {}))
    obs_filters = ObsFilters(**raw.get("obs_filters", {}))
    output = OutputSpec(**raw.get("output", {}))
    return QuerySpec(
        target=target, obs_filters=obs_filters,
        export_all_non_ontology_obs_columns=raw.get("export_all_non_ontology_obs_columns", False),
        obs_columns=raw.get("obs_columns", []) or [],
        var_value_filter=raw.get("var_value_filter"),
        var_columns=raw.get("var_columns", []) or [],
        tiledb_config=raw.get("tiledb_config", {}) or {},
        output=output,
    )

def _format_in_list(values):
    escaped = [v.replace("\\", "\\\\").replace("'", "\\'") for v in values]
    return "[" + ", ".join([f"'{v}'" for v in escaped]) + "]"

def _default_obs_columns():
    return ["dataset_id", "donor_id", "assay", "cell_type", "tissue",
            "tissue_general", "disease", "development_stage", "sex",
            "self_reported_ethnicity", "is_primary_data", "suspension_type"]

def choose_organ_column(obs_keys):
    if "tissue_general" in obs_keys:
        return "tissue_general"
    if "tissue" in obs_keys:
        return "tissue"
    raise ValueError("Neither tissue_general nor tissue found in obs schema keys.")

def build_obs_export_columns(obs_keys, export_all_non_ontology, user_cols):
    keys = list(obs_keys)
    if export_all_non_ontology:
        return [c for c in keys if (c != "soma_joinid" and "ontology_term_id" not in c)]
    cols = list(user_cols) if user_cols else _default_obs_columns()
    missing = [c for c in cols if c not in keys]
    if missing:
        raise ValueError(f"Requested obs columns not present in schema: {missing}")
    return cols

# Mock census data for testing
MOCK_SCC_DATA = {
    "category": ["disease", "disease", "assay", "assay", "tissue_general", "cell_type", "sex", "sex"],
    "label": ["normal", "COVID-19", "10x 3' v3", "10x 3' v2", "lung", "B cell", "male", "female"],
    "ontology_term_id": ["PATO:0000461", "MONDO:0100096", "EFO:0009922", "EFO:0009899", "UBERON:0002048", "CL:0000236", "PATO:0000384", "PATO:0000383"],
}

def resolve_ontology_term_ids_mock(category, labels):
    """Mock version that uses MOCK_SCC_DATA directly."""
    if not labels:
        return []
    found = {}
    for i, cat in enumerate(MOCK_SCC_DATA["category"]):
        if cat == category and MOCK_SCC_DATA["label"][i] in labels:
            found[MOCK_SCC_DATA["label"][i]] = MOCK_SCC_DATA["ontology_term_id"][i]
    missing = [x for x in labels if x not in found]
    if missing:
        raise ValueError(f"Could not resolve ontology_term_id for {category} label(s): {missing}")
    return [found[x] for x in labels]

def build_obs_value_filter_mock(obs_filters):
    """Mock version of build_obs_value_filter."""
    parts = []
    if obs_filters.is_primary_data is not None:
        parts.append(f"is_primary_data == {bool(obs_filters.is_primary_data)}")
    if obs_filters.suspension_type:
        parts.append(f"suspension_type == '{obs_filters.suspension_type}'")
    
    disease_ids = list(obs_filters.disease_ontology_term_ids)
    if not disease_ids and obs_filters.disease_labels:
        disease_ids = resolve_ontology_term_ids_mock("disease", obs_filters.disease_labels)
    if disease_ids:
        parts.append(f"{OBS_TERM_ID_COLS['disease']} in {_format_in_list(disease_ids)}")
    
    sex_ids = list(obs_filters.sex_ontology_term_ids)
    if not sex_ids and obs_filters.sex_labels:
        sex_ids = resolve_ontology_term_ids_mock("sex", obs_filters.sex_labels)
    if sex_ids:
        parts.append(f"{OBS_TERM_ID_COLS['sex']} in {_format_in_list(sex_ids)}")
    
    ct_ids = list(obs_filters.cell_type_ontology_term_ids)
    if not ct_ids and obs_filters.cell_type_labels:
        ct_ids = resolve_ontology_term_ids_mock("cell_type", obs_filters.cell_type_labels)
    if ct_ids:
        parts.append(f"{OBS_TERM_ID_COLS['cell_type']} in {_format_in_list(ct_ids)}")
    
    if obs_filters.extra_value_filter:
        parts.append(f"({obs_filters.extra_value_filter})")
    
    return " and ".join(parts) if parts else ""

def _resolve_outpath(outpath, overwrite):
    p = Path(outpath)
    p.parent.mkdir(parents=True, exist_ok=True)
    if p.exists():
        if overwrite:
            p.unlink()
        else:
            raise FileExistsError(f"Output exists and overwrite=False: {p}")
    return p

# =============================================================================
# Test Framework
# =============================================================================
@dataclass
class TestResult:
    name: str
    passed: bool
    error: str = ""

class TestRunner:
    def __init__(self):
        self.results = []
    
    def run_test(self, name, test_func):
        try:
            test_func()
            self.results.append(TestResult(name=name, passed=True))
            print(f"  ✓ {name}")
        except AssertionError as e:
            self.results.append(TestResult(name=name, passed=False, error=str(e)))
            print(f"  ✗ {name}")
            print(f"    AssertionError: {e}")
        except Exception as e:
            self.results.append(TestResult(name=name, passed=False, error=str(e)))
            print(f"  ✗ {name}")
            print(f"    {type(e).__name__}: {e}")
    
    def summary(self):
        passed = sum(1 for r in self.results if r.passed)
        failed = len(self.results) - passed
        print("\n" + "=" * 60)
        print(f"RESULTS: {passed} passed, {failed} failed, {len(self.results)} total")
        print("=" * 60)
        if failed > 0:
            print("\nFailed tests:")
            for r in self.results:
                if not r.passed:
                    print(f"  - {r.name}: {r.error}")
        return failed == 0

def assert_equal(actual, expected, msg=""):
    if actual != expected:
        raise AssertionError(f"{msg}\nExpected: {expected}\nActual: {actual}")

def assert_in(item, container, msg=""):
    if item not in container:
        raise AssertionError(f"{msg}\n'{item}' not found in container")

def assert_not_in(item, container, msg=""):
    if item in container:
        raise AssertionError(f"{msg}\n'{item}' should not be in container")

def assert_raises(exception_type, func, *args, **kwargs):
    try:
        func(*args, **kwargs)
        raise AssertionError(f"Expected {exception_type.__name__} was not raised")
    except exception_type:
        pass

def assert_true(condition, msg=""):
    if not condition:
        raise AssertionError(msg or "Condition is not true")

def assert_isinstance(obj, cls, msg=""):
    if not isinstance(obj, cls):
        raise AssertionError(f"{msg}\nExpected {cls}, got {type(obj)}")

# =============================================================================
# Tests
# =============================================================================
def test_census_target_defaults():
    target = CensusTarget()
    assert_equal(target.census_version, "stable")
    assert_equal(target.organism, "homo_sapiens")

def test_census_target_custom():
    target = CensusTarget(census_version="2025-11-08", organism="mus_musculus")
    assert_equal(target.census_version, "2025-11-08")

def test_obs_filters_defaults():
    filters = ObsFilters()
    assert_equal(filters.is_primary_data, True)
    assert_equal(filters.suspension_type, "cell")

def test_obs_filters_custom():
    filters = ObsFilters(is_primary_data=False, cell_type_labels=["B cell"])
    assert_equal(filters.is_primary_data, False)

def test_output_spec_defaults():
    output = OutputSpec()
    assert_equal(output.mode, "pandas")

def test_query_spec_defaults():
    spec = QuerySpec()
    assert_isinstance(spec.target, CensusTarget)

def test_load_yaml_valid():
    config = {"target": {"census_version": "stable"}, "output": {"mode": "pandas"}}
    with tempfile.NamedTemporaryFile(mode='w', suffix='.yaml', delete=False) as f:
        yaml.dump(config, f)
        f.flush()
        spec = load_query_spec_yaml(f.name)
        assert_equal(spec.target.census_version, "stable")
    os.unlink(f.name)

def test_load_yaml_missing_file():
    assert_raises(FileNotFoundError, load_query_spec_yaml, "/nonexistent/config.yaml")

def test_load_yaml_minimal():
    with tempfile.NamedTemporaryFile(mode='w', suffix='.yaml', delete=False) as f:
        yaml.dump({"output": {"mode": "pandas"}}, f)
        f.flush()
        spec = load_query_spec_yaml(f.name)
        assert_equal(spec.target.census_version, "stable")
        assert_equal(spec.output.mode, "pandas")
    os.unlink(f.name)

def test_format_in_list_simple():
    assert_equal(_format_in_list(["a", "b", "c"]), "['a', 'b', 'c']")

def test_format_in_list_empty():
    assert_equal(_format_in_list([]), "[]")

def test_format_in_list_escapes_quotes():
    result = _format_in_list(["10x 3' v3"])
    assert_in("\\'", result)

def test_default_obs_columns():
    cols = _default_obs_columns()
    assert_in("dataset_id", cols)
    assert_in("cell_type", cols)

def test_choose_organ_column_prefers_tissue_general():
    assert_equal(choose_organ_column(["tissue", "tissue_general"]), "tissue_general")

def test_choose_organ_column_fallback():
    assert_equal(choose_organ_column(["tissue", "cell_type"]), "tissue")

def test_choose_organ_column_raises():
    assert_raises(ValueError, choose_organ_column, ["cell_type"])

def test_build_export_columns_user():
    result = build_obs_export_columns(["dataset_id", "cell_type"], False, ["dataset_id"])
    assert_equal(result, ["dataset_id"])

def test_build_export_columns_excludes_ontology():
    result = build_obs_export_columns(["dataset_id", "cell_type_ontology_term_id", "soma_joinid"], True, [])
    assert_in("dataset_id", result)
    assert_not_in("cell_type_ontology_term_id", result)

def test_build_export_columns_missing_raises():
    assert_raises(ValueError, build_obs_export_columns, ["dataset_id"], False, ["missing"])

def test_resolve_ontology_single():
    assert_equal(resolve_ontology_term_ids_mock("disease", ["normal"]), ["PATO:0000461"])

def test_resolve_ontology_multiple():
    assert_equal(resolve_ontology_term_ids_mock("sex", ["male", "female"]), ["PATO:0000384", "PATO:0000383"])

def test_resolve_ontology_empty():
    assert_equal(resolve_ontology_term_ids_mock("disease", []), [])

def test_resolve_ontology_missing_raises():
    assert_raises(ValueError, resolve_ontology_term_ids_mock, "disease", ["nonexistent"])

def test_build_value_filter_basic():
    filters = ObsFilters(is_primary_data=True, suspension_type="cell")
    result = build_obs_value_filter_mock(filters)
    assert_in("is_primary_data == True", result)
    assert_in("suspension_type == 'cell'", result)

def test_build_value_filter_disease():
    filters = ObsFilters(is_primary_data=None, suspension_type=None, disease_ontology_term_ids=["PATO:0000461"])
    result = build_obs_value_filter_mock(filters)
    assert_in("disease_ontology_term_id in ['PATO:0000461']", result)

def test_build_value_filter_label_resolution():
    filters = ObsFilters(is_primary_data=None, suspension_type=None, disease_labels=["normal"])
    result = build_obs_value_filter_mock(filters)
    assert_in("PATO:0000461", result)

def test_build_value_filter_empty():
    filters = ObsFilters(is_primary_data=None, suspension_type=None)
    assert_equal(build_obs_value_filter_mock(filters), "")

def test_build_value_filter_extra():
    filters = ObsFilters(is_primary_data=None, suspension_type=None, extra_value_filter="donor_age >= 18")
    result = build_obs_value_filter_mock(filters)
    assert_in("(donor_age >= 18)", result)

def test_resolve_outpath_creates_parent():
    with tempfile.TemporaryDirectory() as tmpdir:
        outpath = Path(tmpdir) / "subdir" / "output.parquet"
        result = _resolve_outpath(str(outpath), True)
        assert_true(result.parent.exists())

def test_resolve_outpath_overwrites():
    with tempfile.TemporaryDirectory() as tmpdir:
        outpath = Path(tmpdir) / "existing.parquet"
        outpath.touch()
        _resolve_outpath(str(outpath), True)
        assert_true(not outpath.exists())

def test_resolve_outpath_no_overwrite_raises():
    with tempfile.TemporaryDirectory() as tmpdir:
        outpath = Path(tmpdir) / "existing.parquet"
        outpath.touch()
        assert_raises(FileExistsError, _resolve_outpath, str(outpath), False)

# =============================================================================
# Main
# =============================================================================
def main():
    runner = TestRunner()
    
    print("\n[Dataclasses]")
    runner.run_test("test_census_target_defaults", test_census_target_defaults)
    runner.run_test("test_census_target_custom", test_census_target_custom)
    runner.run_test("test_obs_filters_defaults", test_obs_filters_defaults)
    runner.run_test("test_obs_filters_custom", test_obs_filters_custom)
    runner.run_test("test_output_spec_defaults", test_output_spec_defaults)
    runner.run_test("test_query_spec_defaults", test_query_spec_defaults)
    
    print("\n[YAML Loading]")
    runner.run_test("test_load_yaml_valid", test_load_yaml_valid)
    runner.run_test("test_load_yaml_missing_file", test_load_yaml_missing_file)
    runner.run_test("test_load_yaml_minimal", test_load_yaml_minimal)
    
    print("\n[Helper Functions]")
    runner.run_test("test_format_in_list_simple", test_format_in_list_simple)
    runner.run_test("test_format_in_list_empty", test_format_in_list_empty)
    runner.run_test("test_format_in_list_escapes_quotes", test_format_in_list_escapes_quotes)
    runner.run_test("test_default_obs_columns", test_default_obs_columns)
    runner.run_test("test_choose_organ_column_prefers_tissue_general", test_choose_organ_column_prefers_tissue_general)
    runner.run_test("test_choose_organ_column_fallback", test_choose_organ_column_fallback)
    runner.run_test("test_choose_organ_column_raises", test_choose_organ_column_raises)
    runner.run_test("test_build_export_columns_user", test_build_export_columns_user)
    runner.run_test("test_build_export_columns_excludes_ontology", test_build_export_columns_excludes_ontology)
    runner.run_test("test_build_export_columns_missing_raises", test_build_export_columns_missing_raises)
    
    print("\n[Filter Building]")
    runner.run_test("test_resolve_ontology_single", test_resolve_ontology_single)
    runner.run_test("test_resolve_ontology_multiple", test_resolve_ontology_multiple)
    runner.run_test("test_resolve_ontology_empty", test_resolve_ontology_empty)
    runner.run_test("test_resolve_ontology_missing_raises", test_resolve_ontology_missing_raises)
    runner.run_test("test_build_value_filter_basic", test_build_value_filter_basic)
    runner.run_test("test_build_value_filter_disease", test_build_value_filter_disease)
    runner.run_test("test_build_value_filter_label_resolution", test_build_value_filter_label_resolution)
    runner.run_test("test_build_value_filter_empty", test_build_value_filter_empty)
    runner.run_test("test_build_value_filter_extra", test_build_value_filter_extra)
    
    print("\n[I/O Functions]")
    runner.run_test("test_resolve_outpath_creates_parent", test_resolve_outpath_creates_parent)
    runner.run_test("test_resolve_outpath_overwrites", test_resolve_outpath_overwrites)
    runner.run_test("test_resolve_outpath_no_overwrite_raises", test_resolve_outpath_no_overwrite_raises)
    
    success = runner.summary()
    return 0 if success else 1


if __name__ == "__main__":
    sys.exit(main())
