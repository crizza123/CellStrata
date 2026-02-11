# census_query Tutorial

A practical guide to querying the CELLxGENE Census with CellStrata's
`census_query` package.

## Prerequisites

Install the required dependencies:

```bash
pip install -r requirements.txt
```

This installs `cellxgene-census`, `pandas`, `pyarrow`, `anndata`, `scanpy`,
`pyyaml`, and other dependencies.

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

`census_query` supports five output modes, each suited for different
downstream workflows. Set the mode in your YAML config or in `OutputSpec`.

### 3a. `pandas` — DataFrame for exploration

Best for interactive analysis, quick inspection, and small-to-medium queries.

```yaml
output:
  mode: pandas
```

```python
df = run_query(spec)               # returns pd.DataFrame
df.groupby("cell_type").size()     # cell type counts
df["donor_id"].nunique()           # number of unique donors
```

### 3b. `arrow` — Columnar format for large datasets

More memory-efficient than pandas. Good for filtering and transformations
before converting to pandas.

```yaml
output:
  mode: arrow
```

```python
table = run_query(spec)            # returns pyarrow.Table
print(f"{table.num_rows:,} rows, {table.nbytes / 1e6:.1f} MB")
# Convert a subset to pandas when ready:
subset = table.filter(table.column("sex") == "female").to_pandas()
```

### 3c. `anndata` — Expression matrix for single-cell analysis

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

### 3d. `parquet` — Stream large results to disk

For very large queries that don't fit in memory. Writes results in batches.

```yaml
output:
  mode: parquet
  outpath: output/healthy_cells.parquet
  parquet_compression: zstd       # or snappy, gzip, none
  overwrite: true
```

```python
outpath = run_query(spec)          # returns Path to the file
print(f"Written to: {outpath}")

# Read back later:
import pandas as pd
df = pd.read_parquet(outpath)
```

### 3e. `dataset_list` — Find matching CELLxGENE datasets

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

When `outpath` is set, the IDs are also written to a text file (one per line)
that can be piped into other tools:

```bash
# Count matching datasets
wc -l output/matching_datasets.txt

# Use with another script
while read ds_id; do
    echo "Processing dataset: $ds_id"
done < output/matching_datasets.txt
```

---

## 4. Filter Reference

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

## 5. Practical Examples

### Example A: Compare mast cells across tissues

```python
from census_query import QuerySpec, ObsFilters, OutputSpec, run_query

spec = QuerySpec(
    obs_filters=ObsFilters(
        is_primary_data=True,
        disease_ontology_term_ids=["PATO:0000461"],
        cell_type_labels=["mast cell"],
    ),
    output=OutputSpec(mode="pandas"),
)

df = run_query(spec)

# How many mast cells per tissue?
print(df.groupby("tissue_general")["donor_id"].agg(
    cells="size",
    donors="nunique",
).sort_values("cells", ascending=False))
```

### Example B: Export healthy lung data to Parquet for later analysis

```yaml
# lung_export.yaml
target:
  census_version: stable
  organism: homo_sapiens

obs_filters:
  is_primary_data: true
  suspension_type: cell
  disease_ontology_term_ids: ["PATO:0000461"]
  tissue_general_labels: ["lung"]
  assay_labels: ["10x 3' v3", "10x 3' v2"]

output:
  mode: parquet
  outpath: output/healthy_lung_metadata.parquet
  overwrite: true
```

```bash
python -m census_query --config lung_export.yaml
# Output: Wrote Parquet: 1,234,567 rows in 25 batches to output/healthy_lung_metadata.parquet
```

### Example C: Find datasets containing specific cell types

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
    print(f"  https://cellxgene.cziscience.com/e/{ds_id}.cxg/")
```

### Example D: Download expression data for a scanpy pipeline

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

## 6. S3 Timeout Tuning

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

## 7. Package Structure

After the recent refactor, `census_query` is organized as a sub-package:

```
census_query/
  __init__.py       # facade — re-exports all public symbols
  __main__.py       # CLI entry point (python -m census_query)
  _schema.py        # dataclasses, constants, YAML loading
  _filters.py       # filter building, ontology resolution
  _io.py            # streaming, Parquet writing, TileDB context
  _runner.py        # run_query() orchestrator, CLI main()
```

You can import from the top-level package for convenience:

```python
from census_query import run_query, QuerySpec, ObsFilters
```

Or import directly from sub-modules when you only need part of the
functionality:

```python
from census_query._schema import load_query_spec_yaml    # lightweight, no pandas
from census_query._filters import resolve_ontology_term_ids  # filters only
```
