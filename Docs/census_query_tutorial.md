# census_query Tutorial

A practical guide to querying the CELLxGENE Census and visualizing cell
metadata with CellStrata's `census_query` package.

## Prerequisites

Install the required dependencies:

```bash
pip install -r requirements.txt
```

This installs `cellxgene-census`, `pandas`, `anndata`, `scanpy`,
`matplotlib`, `seaborn`, `pyyaml`, and other dependencies.

---

## 1. Quick Start — Run a Query from a YAML Config

The fastest way to use `census_query` is with the CLI and a YAML config file.

### Step 1: Write a config file

Create `my_query.yaml`:

```yaml
target:
  census_version: stable
  organism: homo_sapiens
  measurement: RNA

obs_filters:
  is_primary_data: true
  suspension_type: cell
  disease_ontology_term_ids: ["PATO:0000461"]   # healthy/normal
  assay_labels: ["10x 3' v3"]
  tissue_general_labels: ["blood"]

output:
  mode: pandas
```

### Step 2: Run the query

```bash
python -m census_query --config my_query.yaml
```

Output:

```
Query returned DataFrame with 245,893 rows
Columns: ['dataset_id', 'donor_id', 'assay', 'cell_type', 'tissue',
          'tissue_general', 'disease', 'development_stage', 'sex',
          'self_reported_ethnicity', 'is_primary_data', 'suspension_type']

First 5 rows:
                           dataset_id    donor_id      assay cell_type ...
0  a1b2c3d4-e5f6-...  donor_001  10x 3' v3    T cell ...
```

Each row is one **cell** from the Census. The CLI prints a summary and the
first 5 rows so you can verify your filters are working as expected.

---

## 2. Using census_query as a Python Module

For more control, import `census_query` directly in your scripts or notebooks.

### 2a. Load a YAML config and run a query

```python
from census_query import load_query_spec_yaml, run_query

spec = load_query_spec_yaml("config/census_query.yaml")
df = run_query(spec)

print(f"{len(df):,} cells returned")
print(df["cell_type"].value_counts().head(10))
```

### 2b. Build a query entirely in Python (no YAML)

You can skip the YAML file and construct the query programmatically:

```python
from census_query import (
    CensusTarget,
    ObsFilters,
    OutputSpec,
    QuerySpec,
    run_query,
)

spec = QuerySpec(
    target=CensusTarget(
        census_version="stable",
        organism="homo_sapiens",
    ),
    obs_filters=ObsFilters(
        is_primary_data=True,
        suspension_type="cell",
        disease_ontology_term_ids=["PATO:0000461"],
        assay_labels=["10x 3' v3"],
        tissue_general_labels=["lung"],
        cell_type_labels=["mast cell"],
    ),
    output=OutputSpec(mode="pandas"),
)

df = run_query(spec)
print(f"Found {len(df):,} mast cells in healthy lung tissue")
print(df[["dataset_id", "donor_id", "cell_type", "tissue"]].head())
```

This is useful when you want to parameterize queries in a loop or integrate
`census_query` into a larger pipeline.

---

## 3. Output Modes

`census_query` supports three output modes. Set the mode in your YAML config
or in `OutputSpec`.

### 3a. `pandas` — DataFrame for exploration and visualization

Best for interactive analysis, quick inspection, and metadata visualization.

```yaml
output:
  mode: pandas
```

```python
df = run_query(spec)               # returns pd.DataFrame
df.groupby("cell_type").size()     # cell type counts
df["donor_id"].nunique()           # number of unique donors
```

### 3b. `anndata` — Expression matrix for single-cell analysis

This is the only mode that downloads **gene expression data**. Returns an
`AnnData` object ready for use with scanpy.

```yaml
output:
  mode: anndata

# Optional: filter to specific genes
var_value_filter: "feature_name in ['TPSAB1', 'TPSB2', 'CMA1', 'KIT']"
```

```python
import scanpy as sc

adata = run_query(spec)
print(f"{adata.n_obs} cells x {adata.n_vars} genes")

# Standard scanpy workflow
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata)
sc.tl.pca(adata)
sc.tl.umap(adata)
sc.pl.umap(adata, color="cell_type")
```

### 3c. `dataset_list` — Find matching CELLxGENE datasets

Returns a sorted list of unique dataset IDs whose cells match your filters.
Useful when you want to identify which datasets to download from the
CELLxGENE portal without pulling all the metadata.

```yaml
output:
  mode: dataset_list
  outpath: output/matching_datasets.txt   # optional
```

```python
dataset_ids = run_query(spec)      # returns List[str]
print(f"{len(dataset_ids)} datasets contain matching cells:")
for ds_id in dataset_ids:
    print(f"  {ds_id}")
```

When `outpath` is set, the IDs are also written to a text file (one per line).

---

## 4. Visualizing Metadata

After querying cell metadata in `pandas` mode, use the `_visualize` module
to explore distributions with publication-ready plots.

### 4a. Individual plots

Every plot function accepts a `pd.DataFrame` and returns a `matplotlib.Axes`,
so you can compose and customize figures easily.

```python
from census_query import (
    run_query,
    load_query_spec_yaml,
    plot_cell_type_counts,
    plot_tissue_counts,
    plot_sex_distribution,
    plot_assay_counts,
    plot_disease_counts,
    plot_dataset_contribution,
    plot_development_stage_counts,
    plot_donors_per_dataset,
    plot_cell_type_by_tissue,
)
import matplotlib.pyplot as plt

# 1. Run a query
spec = load_query_spec_yaml("config/census_query.yaml")
df = run_query(spec)

# 2. Plot individual distributions
plot_cell_type_counts(df, top_n=15)
plt.tight_layout()
plt.show()

plot_tissue_counts(df, column="tissue_general")
plt.tight_layout()
plt.show()

plot_sex_distribution(df)
plt.tight_layout()
plt.show()
```

### 4b. Composing a custom multi-panel figure

Because every function accepts an `ax` parameter, you can build your own
layouts:

```python
fig, axes = plt.subplots(1, 3, figsize=(18, 5))

plot_cell_type_counts(df, top_n=10, ax=axes[0])
plot_tissue_counts(df, ax=axes[1])
plot_sex_distribution(df, ax=axes[2])

fig.suptitle("My custom metadata overview", fontsize=14)
fig.tight_layout()
fig.savefig("custom_overview.png", dpi=150)
```

### 4c. Cell-type-by-tissue heatmap

The heatmap shows how cell types distribute across tissues. Colour encodes
`log10(count + 1)` so rare populations remain visible; cell annotations
show raw counts.

```python
plot_cell_type_by_tissue(df, top_cell_types=12, top_tissues=8)
plt.tight_layout()
plt.show()
```

### 4d. One-line summary dashboard

`plot_metadata_summary` creates a 3x2 grid with six panels:

1. Cell type counts
2. Tissue counts
3. Sex distribution
4. Assay distribution
5. Dataset contribution
6. Cell-type-by-tissue heatmap

```python
from census_query import plot_metadata_summary

fig = plot_metadata_summary(df, top_n=15, save_path="figures/summary.png")
plt.show()
```

The dashboard title automatically includes the total cell count, unique
donors, and unique datasets.

---

## 5. Available Plot Functions

| Function | What it shows | Key parameters |
|---|---|---|
| `plot_cell_type_counts(df)` | Horizontal bar chart of cell-type frequencies | `top_n`, `ax` |
| `plot_tissue_counts(df)` | Horizontal bar chart of tissue frequencies | `column`, `top_n`, `ax` |
| `plot_sex_distribution(df)` | Bar chart of sex distribution | `ax` |
| `plot_assay_counts(df)` | Bar chart of sequencing assay frequencies | `ax` |
| `plot_disease_counts(df)` | Horizontal bar chart of disease annotations | `top_n`, `ax` |
| `plot_dataset_contribution(df)` | Cells per dataset (horizontal bar) | `top_n`, `ax` |
| `plot_development_stage_counts(df)` | Bar chart of development stages | `top_n`, `ax` |
| `plot_donors_per_dataset(df)` | Unique donors per dataset (horizontal bar) | `top_n`, `ax` |
| `plot_cell_type_by_tissue(df)` | Heatmap: cell types vs tissues | `top_cell_types`, `top_tissues`, `tissue_column`, `ax` |
| `plot_metadata_summary(df)` | 3x2 dashboard combining six panels | `top_n`, `tissue_column`, `save_path`, `dpi` |

All individual functions return `plt.Axes`; `plot_metadata_summary` returns
`plt.Figure`.

---

## 6. Filter Reference

All filters go under `obs_filters` in the YAML config. You can use
human-readable **labels** (automatically resolved to ontology term IDs) or
provide **ontology term IDs** directly for stability.

### Labels vs. ontology term IDs

```yaml
obs_filters:
  # Option A: human-readable labels (resolved at query time)
  disease_labels: ["normal", "COVID-19"]
  assay_labels: ["10x 3' v3"]
  tissue_general_labels: ["lung", "blood"]
  cell_type_labels: ["mast cell", "T cell"]
  sex_labels: ["female"]

  # Option B: ontology term IDs (stable across Census versions)
  disease_ontology_term_ids: ["PATO:0000461"]
  assay_ontology_term_ids: ["EFO:0009922"]
  tissue_general_ontology_term_ids: ["UBERON:0002048"]
  cell_type_ontology_term_ids: ["CL:0000097"]
  sex_ontology_term_ids: ["PATO:0000383"]
```

If both labels and ontology term IDs are provided for the same category,
the ontology term IDs take precedence.

### Boolean and string filters

```yaml
obs_filters:
  is_primary_data: true          # exclude duplicate/re-annotated cells
  suspension_type: cell          # "cell" or "nucleus"
```

### Custom filter expression

For filters not covered by the built-in fields, use `extra_value_filter`
with raw SOMA query syntax:

```yaml
obs_filters:
  extra_value_filter: "development_stage != 'unknown'"
```

### Available filter categories

| Category | Label field | Ontology ID field | Example label |
|----------|------------|-------------------|--------------|
| Disease | `disease_labels` | `disease_ontology_term_ids` | `"normal"`, `"COVID-19"` |
| Sex | `sex_labels` | `sex_ontology_term_ids` | `"male"`, `"female"` |
| Assay | `assay_labels` | `assay_ontology_term_ids` | `"10x 3' v3"`, `"10x multiome"` |
| Tissue (coarse) | `tissue_general_labels` | `tissue_general_ontology_term_ids` | `"lung"`, `"blood"`, `"brain"` |
| Tissue (fine) | `tissue_labels` | `tissue_ontology_term_ids` | `"lung parenchyma"` |
| Cell type | `cell_type_labels` | `cell_type_ontology_term_ids` | `"mast cell"`, `"B cell"` |
| Dev. stage | `development_stage_labels` | `development_stage_ontology_term_ids` | `"adult"` |

---

## 7. Practical Examples

### Example A: Compare mast cells across tissues and visualize

```python
from census_query import (
    QuerySpec, ObsFilters, OutputSpec, run_query,
    plot_cell_type_by_tissue, plot_metadata_summary,
)

spec = QuerySpec(
    obs_filters=ObsFilters(
        is_primary_data=True,
        disease_ontology_term_ids=["PATO:0000461"],
        cell_type_labels=["mast cell"],
    ),
    output=OutputSpec(mode="pandas"),
)

df = run_query(spec)

# Quick tabular summary
print(df.groupby("tissue_general")["donor_id"].agg(
    cells="size",
    donors="nunique",
).sort_values("cells", ascending=False))

# Visual summary
fig = plot_metadata_summary(df, save_path="mast_cell_summary.png")
```

### Example B: Find datasets containing specific cell types

```python
from census_query import QuerySpec, ObsFilters, OutputSpec, run_query

spec = QuerySpec(
    obs_filters=ObsFilters(
        is_primary_data=True,
        disease_ontology_term_ids=["PATO:0000461"],
        cell_type_labels=["mast cell"],
        tissue_general_labels=["lung"],
    ),
    output=OutputSpec(
        mode="dataset_list",
        outpath="output/mast_cell_lung_datasets.txt",
    ),
)

dataset_ids = run_query(spec)
print(f"{len(dataset_ids)} datasets have mast cells in healthy lung:")
for ds_id in dataset_ids:
    print(f"  {ds_id}")
```

### Example C: Download expression data for a scanpy pipeline

```python
from census_query import QuerySpec, ObsFilters, OutputSpec, run_query
import scanpy as sc

spec = QuerySpec(
    obs_filters=ObsFilters(
        is_primary_data=True,
        disease_ontology_term_ids=["PATO:0000461"],
        cell_type_labels=["mast cell"],
        tissue_general_labels=["lung"],
        assay_labels=["10x 3' v3"],
    ),
    output=OutputSpec(mode="anndata"),
    var_value_filter="feature_name in ['TPSAB1', 'TPSB2', 'CMA1', 'KIT', 'FCER1A']",
)

adata = run_query(spec)

# Quick QC
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

# Visualize marker expression across donors
sc.pl.dotplot(adata, var_names=["TPSAB1", "TPSB2", "CMA1", "KIT"], groupby="donor_id")
```

### Example D: Custom two-panel comparison figure

```python
import matplotlib.pyplot as plt
from census_query import (
    QuerySpec, ObsFilters, OutputSpec, run_query,
    plot_cell_type_counts, plot_donors_per_dataset,
)

spec = QuerySpec(
    obs_filters=ObsFilters(
        is_primary_data=True,
        disease_ontology_term_ids=["PATO:0000461"],
        tissue_general_labels=["lung"],
    ),
    output=OutputSpec(mode="pandas"),
)

df = run_query(spec)

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 6))
plot_cell_type_counts(df, top_n=10, ax=ax1)
plot_donors_per_dataset(df, top_n=10, ax=ax2)
fig.suptitle("Healthy lung: cell types and donor coverage")
fig.tight_layout()
fig.savefig("lung_overview.png", dpi=150)
```

### Example E: Using the filter module independently

Because `census_query` is split into independent sub-modules, you can use
the filter builder on its own — for example, to inspect what filter string
would be generated without actually running a query:

```python
import cellxgene_census as cxc
from census_query._filters import build_obs_value_filter
from census_query import ObsFilters

filters = ObsFilters(
    is_primary_data=True,
    disease_ontology_term_ids=["PATO:0000461"],
    assay_labels=["10x 3' v3"],
    cell_type_labels=["mast cell"],
)

with cxc.open_soma() as census:
    filter_string = build_obs_value_filter(census, filters)
    print(filter_string)
    # is_primary_data == True and suspension_type == 'cell'
    #   and disease_ontology_term_id in ['PATO:0000461']
    #   and assay_ontology_term_id in ['EFO:0009922']
    #   and cell_type_ontology_term_id in ['CL:0000097']
```

---

## 8. S3 Timeout Tuning

Large Census queries can sometimes hit S3 timeouts. Add `tiledb_config` to
your YAML to mitigate this:

```yaml
tiledb_config:
  vfs.s3.connect_timeout_ms: 60000      # 60 seconds
  vfs.s3.request_timeout_ms: 600000     # 10 minutes
  vfs.s3.max_parallel_ops: 2            # reduce parallelism
```

This is optional — only needed if you encounter timeout errors on large
queries.

---

## 9. Package Structure

```
census_query/
  __init__.py       # facade — re-exports all public symbols
  __main__.py       # CLI entry point (python -m census_query)
  _schema.py        # dataclasses, constants, YAML loading
  _filters.py       # filter building, ontology resolution
  _runner.py        # run_query() orchestrator, CLI main()
  _visualize.py     # metadata plotting (matplotlib/seaborn)
```

You can import from the top-level package for convenience:

```python
from census_query import run_query, QuerySpec, plot_metadata_summary
```

Or import directly from sub-modules when you only need part of the
functionality:

```python
from census_query._schema import load_query_spec_yaml    # lightweight, no pandas
from census_query._filters import resolve_ontology_term_ids  # filters only
from census_query._visualize import plot_cell_type_counts  # single plot
```
