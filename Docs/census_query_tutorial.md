# census_query Tutorial

A practical guide to querying the CELLxGENE Census with CellStrata's
`census` package. The running example throughout this tutorial is
**T cells in healthy human brain**.

## Prerequisites

Install the package (editable mode recommended during development):

```bash
pip install -e ".[dev]"
```

This installs `cellxgene-census`, `pandas`, `pyarrow`, `anndata`, `scanpy`,
`matplotlib`, `seaborn`, `pyyaml`, and other dependencies.

Alternatively, use conda:

```bash
conda env create -f envs/environment.yml
conda activate cellstrata
pip install -e ".[dev]"
```

---

## 1. Quick Start — Run a Query from a YAML Config

The fastest way to use `census_query` is with the CLI script and a YAML
config file.

### Step 1: Write a config file

Create `my_query.yaml` (or use the included example at
`CellStrata/census_query.yaml`):

```yaml
target:
  census_version: stable
  organism: homo_sapiens

obs_filters:
  is_primary_data: true
  disease: ["normal"]
  tissue_general: ["brain"]
  cell_type: ["T cell"]
  assay: ["10x 3' v2", "10x 3' v3", "10x multiome"]
  sex: ["male", "female"]

output:
  summary_csv: results.csv
  # celllevel_csv: cell_data.csv
  # write_celllevel: false
```

### Step 2: Run the query

```bash
python CellStrata/scripts/run_census_query.py --config CellStrata/census_query.yaml
```

Or, if you installed the package with `pip install -e .`:

```bash
cellstrata-query --config CellStrata/census_query.yaml
```

Output:

```
VALUE_FILTER: is_primary_data == True and disease == 'normal' and tissue_general == 'brain' and cell_type == 'T cell' and assay in ['10x 3' v2', '10x 3' v3', '10x multiome'] and sex in ['male', 'female']
Got 4,821 cells
Wrote results.csv (12 datasets)
```

The script writes a `results.csv` summarizing the number of cells per
dataset. Each row is one dataset with its cell count.

---

## 2. Using the Census Package as a Python Module

For more control, import the census modules directly in your scripts or
notebooks.

### 2a. Load a YAML config and run a query

```python
from CellStrata.census.config import load_config
from CellStrata.census.query import fetch_obs_df, summarize_by_dataset
from CellStrata.census.outputs import write_summary_csv

cfg = load_config("CellStrata/census_query.yaml")
df, value_filter = fetch_obs_df(cfg)

print(f"{len(df):,} cells returned")
print(df["cell_type"].value_counts().head(10))

summary = summarize_by_dataset(df)
write_summary_csv(summary, cfg.output.summary_csv)
```

### 2b. Build a config entirely in Python (no YAML)

You can skip the YAML file and construct the config programmatically:

```python
from CellStrata.census.config import (
    CensusTarget,
    ObsFilters,
    OutputSpec,
    CensusQueryConfig,
)
from CellStrata.census.query import fetch_obs_df, summarize_by_dataset
from CellStrata.census.outputs import write_summary_csv

cfg = CensusQueryConfig(
    target=CensusTarget(
        census_version="stable",
        organism="homo_sapiens",
    ),
    obs_filters=ObsFilters(
        is_primary_data=True,
        disease=["normal"],
        tissue_general=["brain"],
        cell_type=["T cell"],
    ),
    output=OutputSpec(summary_csv="results.csv"),
)

df, value_filter = fetch_obs_df(cfg)
print(f"Found {len(df):,} T cells in healthy brain tissue")
print(df[["dataset_id", "donor_id", "cell_type", "tissue"]].head())

summary = summarize_by_dataset(df)
write_summary_csv(summary, cfg.output.summary_csv)
```

This is useful when you want to parameterize queries in a loop or integrate
the census module into a larger pipeline.

---

## 3. Output Options

The `OutputSpec` controls what gets written to disk.

### 3a. Summary CSV only (default)

By default only a dataset-level summary is written:

```yaml
output:
  summary_csv: results.csv
```

```python
summary = summarize_by_dataset(df)
write_summary_csv(summary, "results.csv")
```

The CSV contains two columns: `dataset_id` and `cell_count`.

### 3b. Cell-level CSV

Enable `write_celllevel` to also dump every individual cell record:

```yaml
output:
  summary_csv: results.csv
  celllevel_csv: cell_data.csv
  write_celllevel: true
```

```python
from CellStrata.census.outputs import write_celllevel_csv

if cfg.output.write_celllevel and cfg.output.celllevel_csv:
    write_celllevel_csv(df, cfg.output.celllevel_csv)
```

---

## 4. Filter Reference

All filters go under `obs_filters` in the YAML config or in the `ObsFilters`
dataclass.

### Available filter fields

| Field | Type | Description | Example |
|-------|------|-------------|---------|
| `is_primary_data` | `bool` | Exclude duplicate/re-annotated cells | `true` |
| `dataset_id` | `str` | Restrict to a single dataset | `"1b350d0a-..."` |
| `disease` | `list[str]` | Filter by disease label | `["normal"]` |
| `tissue_general` | `list[str]` | Filter by coarse tissue label | `["brain", "lung"]` |
| `cell_type` | `list[str]` | Filter by cell type label | `["T cell", "mast cell"]` |
| `assay` | `list[str]` | Filter by sequencing assay | `["10x 3' v3"]` |
| `sex` | `list[str]` | Filter by sex | `["male", "female"]` |
| `extra_value_filter` | `str` | Raw SOMA query clause appended with `and` | `"development_stage != 'unknown'"` |

### Example: combining filters

```yaml
obs_filters:
  is_primary_data: true
  disease: ["normal"]
  tissue_general: ["lung"]
  cell_type: ["mast cell"]
  assay: ["10x 3' v3"]
  sex: ["female"]
  extra_value_filter: "development_stage != 'unknown'"
```

This generates a value filter string like:

```
is_primary_data == True and disease == 'normal' and tissue_general == 'lung'
  and cell_type == 'mast cell' and assay == '10x 3' v3' and sex == 'female'
  and (development_stage != 'unknown')
```

### Inspecting the generated filter

You can inspect the filter string without running the full query:

```python
from CellStrata.census.config import ObsFilters
from CellStrata.census.filters import build_value_filter

filters = ObsFilters(
    is_primary_data=True,
    disease=["normal"],
    tissue_general=["brain"],
    cell_type=["T cell"],
)

print(build_value_filter(filters))
# is_primary_data == True and disease == 'normal' and tissue_general == 'brain' and cell_type == 'T cell'
```

---

## 5. Configuring Observed Columns

By default the query returns these metadata columns:

```
dataset_id, donor_id, assay, sex, tissue, tissue_general,
cell_type, disease, development_stage
```

Override in the YAML with `obs_cols`:

```yaml
obs_cols:
  - dataset_id
  - donor_id
  - cell_type
  - tissue_general
```

Or in Python:

```python
cfg = CensusQueryConfig(
    obs_cols=["dataset_id", "donor_id", "cell_type", "tissue_general"],
    ...
)
```

The script validates that all requested columns exist in the Census schema
before running the query.

---

## 6. Practical Examples

### Example A: Narrow to a single dataset for rapid prototyping

Use `dataset_id` to restrict to one known dataset. This is much faster for
testing and gives you a reproducible, fixed-size slice.

```python
from CellStrata.census.config import CensusQueryConfig, ObsFilters, OutputSpec
from CellStrata.census.query import fetch_obs_df, summarize_by_dataset

cfg = CensusQueryConfig(
    obs_filters=ObsFilters(
        is_primary_data=True,
        dataset_id="1b350d0a-4535-4879-beb6-1142f3f94947",
    ),
    output=OutputSpec(summary_csv="single_dataset.csv"),
)

df, _ = fetch_obs_df(cfg)
print(f"{len(df):,} cells from single dataset")
print(df["cell_type"].value_counts().head(10))
```

### Example B: Cell-level export for downstream analysis

```yaml
# config/full_export.yaml
target:
  census_version: stable
  organism: homo_sapiens

obs_filters:
  is_primary_data: true
  disease: ["normal"]
  tissue_general: ["lung"]
  cell_type: ["mast cell"]

output:
  summary_csv: output/mast_cell_summary.csv
  celllevel_csv: output/mast_cell_data.csv
  write_celllevel: true
```

```bash
python CellStrata/scripts/run_census_query.py --config config/full_export.yaml
```

### Example C: Loop over multiple tissues

```python
from CellStrata.census.config import CensusQueryConfig, ObsFilters, OutputSpec
from CellStrata.census.query import fetch_obs_df, summarize_by_dataset
from CellStrata.census.outputs import write_summary_csv

for tissue in ["lung", "brain", "blood"]:
    cfg = CensusQueryConfig(
        obs_filters=ObsFilters(
            is_primary_data=True,
            disease=["normal"],
            tissue_general=[tissue],
        ),
        output=OutputSpec(summary_csv=f"output/{tissue}_summary.csv"),
    )
    df, _ = fetch_obs_df(cfg)
    summary = summarize_by_dataset(df)
    write_summary_csv(summary, cfg.output.summary_csv)
    print(f"{tissue}: {len(df):,} cells across {len(summary)} datasets")
```

---

## 7. Running Tests

Run the full test suite (no network access required — Census API is mocked):

```bash
python -m pytest CellStrata/tests/test_census_query.py -v
```

Quick smoke test to verify the package imports and config loading:

```bash
python -c "from CellStrata.census.config import load_config; cfg = load_config('CellStrata/config/census_query.yaml'); print(cfg)"
```

---

## 8. Package Structure

```
CellStrata/
  census/
    __init__.py
    config.py       # dataclasses: CensusTarget, ObsFilters, OutputSpec,
                    #   CensusQueryConfig, load_config()
    filters.py      # build_value_filter() — constructs SOMA filter strings
    query.py        # fetch_obs_df(), summarize_by_dataset(), validate_obs_cols()
    outputs.py      # write_summary_csv(), write_celllevel_csv()
  scripts/
    run_census_query.py   # CLI entry point (also: cellstrata-query)
  tests/
    test_census_query.py  # unit tests (mocked, no network)
  config/
    census_query.yaml     # example config
```
