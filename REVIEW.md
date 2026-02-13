# Code Review: `scripts/run_census_query.py` (branch `refactor_1`)

## Overview

The refactor moved from a monolithic `census_query` sub-package (with `_schema`,
`_filters`, `_io`, `_runner`, `_visualize` modules) to a simpler `census` package
with four small files: `config.py`, `filters.py`, `outputs.py`, `query.py`.
The entry-point script is `scripts/run_census_query.py`.

## Verdict: **Not Functional** - several blocking issues prevent execution

---

## CRITICAL Issues (will crash at runtime)

### 1. Wrong import path in `scripts/run_census_query.py` (lines 5-7)

```python
from cellstrata.census.config import load_config
from cellstrata.census.query import fetch_obs_df, summarize_by_dataset
from cellstrata.census.outputs import write_summary_csv, write_celllevel_csv
```

The package directory is `CellStrata` (capital C, capital S) but imports reference
`cellstrata` (all lowercase). There is also no top-level `CellStrata/__init__.py`
on this branch. This will produce `ModuleNotFoundError` immediately.

### 2. YAML configs are incompatible with the new dataclasses

Both `census_query.yaml` and `config/census_query.yaml` still use the **old** schema
fields (`suspension_type`, `disease_ontology_term_ids`, `tissue_general_labels`,
`cell_type_labels`, `assay_labels`, `sex_labels`, `output.mode`, `output.outpath`,
`tiledb_config`).

The new `ObsFilters` dataclass expects `sex`, `disease`, `tissue_general`,
`cell_type`, `assay` (flat label lists). The new `OutputSpec` expects `summary_csv`,
`celllevel_csv`, `write_celllevel`.

Feeding the existing YAMLs to `load_config()` will raise:
`TypeError: __init__() got an unexpected keyword argument 'suspension_type'`

### 3. `query.py` accesses `human.obs.schema.names` (line 11)

```python
schema_names = set(human.obs.schema.names)
```

The `tiledbsoma` `SOMADataFrame` exposes `.keys()` not `.schema.names` in all
Census SDK versions. This may raise `AttributeError` depending on the installed
version. The old code used `exp.obs.keys()`.

### 4. Test file (`test_census_query.py`) is entirely stale

The test file still imports from the **old** `census_query` package:
```python
from census_query import (CensusTarget, ObsFilters, OutputSpec, QuerySpec, ...)
```

None of these exist in the new `census` package. Every test will fail at
collection with `ImportError`. The test suite provides zero coverage for the
refactored code.

### 5. `pyproject.toml` entry point is broken

```toml
cellstrata-query = "CellStrata.census_query._runner:main"
```

`CellStrata.census_query._runner` no longer exists. The entry point should
reference the new script or the new package structure.

---

## MODERATE Issues

### 6. `load_config()` has no error handling

```python
cfg = yaml.safe_load(Path(path).read_text())
```

If the file doesn't exist, you get a raw `FileNotFoundError` with no helpful
message. If YAML is empty, `cfg` is `None` and `cfg.get(...)` raises
`AttributeError`.

### 7. `query.py:validate_obs_cols` opens the experiment twice

`validate_obs_cols` accesses `census["census_data"][organism]`, then `fetch_obs_df`
accesses it again on line 23. Minor but unnecessary duplication.

### 8. No `__init__.py` at the `CellStrata/` package root

Without this, `import CellStrata` or `from CellStrata.census import ...` won't
work unless the project is installed in editable mode with the right setuptools
configuration. The old branch had one.

### 9. Filter stability regression

The old code resolved human-readable labels to ontology term IDs for stable
queries across Census versions. The refactored `filters.py` passes labels directly
(e.g., `disease == 'normal'`), which is less stable if Census renames labels.

---

## MINOR / Cleanup

- **Stale `__pycache__`** in `census_query/` should be removed.
- **`requirements.txt`** duplicates `pyproject.toml` dependencies.
- Both YAML config files are for the old schema and need updating.
- The empty `census_query/` directory (with only stale `.pyc` files) should be
  removed.

---

## Summary

| Category                | Count |
|-------------------------|-------|
| Blocking runtime errors | 5     |
| Moderate issues         | 4     |
| Minor cleanup           | 4     |

The script **will not run** as-is. The refactor appears incomplete - the new
`census/` package is a clean start but the surrounding infrastructure (YAML
configs, tests, pyproject entry point, script imports) was not updated to match.
