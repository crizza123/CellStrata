# census_query - first function tutorial.

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
python CellStrata/scripts/run_census_query.py --config CellStrata/config/census_query.yaml
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
    census_query.yaml     # example config
```
