# Work Breakdown Structure (WBS)

## Activity 1: CellxGene Census API query for data
**Deliverable:** Local dataset(s) retrieved from CellxGene Census (e.g., downloaded `.h5ad` files and/or an extracted dataset subset saved locally)

### Action 1: Define query inputs
- Sub-action 1.1: Identify dataset scope and identifiers to query
  - Task 1.1.1: Specify which dataset identifiers (or selection criteria) will be used for the query
    - Deliverable: Documented list of dataset identifiers or selection rule
    - Completion criteria: List/rule exists and is recorded in the repo (e.g., in config or docs)
- Sub-action 1.2: Define required metadata/fields to retrieve
  - Task 1.2.1: Specify which metadata fields are needed downstream (e.g., dataset_id, tissue)
    - Deliverable: Documented list of required fields
    - Completion criteria: Required fields listed and saved in documentation/config

### Action 2: Retrieve data
- Sub-action 2.1: Execute Census query
  - Task 2.1.1: Implement and run query to retrieve expression + required metadata for the selected scope
    - Deliverable: Retrieved data stored locally
    - Completion criteria: Data can be loaded from disk without errors

### Action 3: Standardize output format for downstream steps
- Sub-action 3.1: Store retrieved data in consistent structure
  - Task 3.1.1: Save output in chosen intermediate format (e.g., `.h5ad`) with consistent naming
    - Deliverable: Saved file(s) in predictable directory structure
    - Completion criteria: Files exist and follow naming/location conventions


## Activity 2: Quality Control
**Deliverable:** QC-filtered dataset stored in the project’s working format (e.g., `.h5ad` with filtered cells/genes)

### Action 1: Compute QC metrics
- Sub-action 1.1: Compute per-cell QC metrics
  - Task 1.1.1: Calculate QC metrics required for filtering
    - Deliverable: Dataset with QC metrics added
    - Completion criteria: QC metric fields exist in the data object
- Sub-action 1.2: Compute per-gene QC metrics (if needed)
  - Task 1.2.1: Calculate gene-level summary metrics used for filtering
    - Deliverable: Dataset with gene-level metrics added
    - Completion criteria: Gene metric fields exist in the data object

### Action 2: Apply QC filters
- Sub-action 2.1: Filter low-quality cells
  - Task 2.1.1: Apply defined filtering rules to remove cells failing QC thresholds
    - Deliverable: Filtered dataset
    - Completion criteria: Filtered dataset size reflects applied rules and loads without errors
- Sub-action 2.2: Filter low-information genes
  - Task 2.2.1: Apply defined rules to remove genes not meeting minimum criteria
    - Deliverable: Filtered dataset
    - Completion criteria: Filtered dataset includes only retained genes and loads without errors

### Action 3: Save QC output
- Sub-action 3.1: Write QC dataset to disk
  - Task 3.1.1: Save QC-filtered dataset
    - Deliverable: QC output file(s)
    - Completion criteria: Output file(s) exist and can be reloaded


## Activity 3: Normalization
**Deliverable:** Normalized dataset saved in project working format

### Action 1: Normalize expression values
- Sub-action 1.1: Apply normalization method
  - Task 1.1.1: Normalize expression values using the selected method
    - Deliverable: Dataset with normalized expression stored
    - Completion criteria: Normalized values are present in the expected location/field

### Action 2: Save normalized output
- Sub-action 2.1: Write normalized dataset to disk
  - Task 2.1.1: Save normalized dataset
    - Deliverable: Normalized output file(s)
    - Completion criteria: Output file(s) exist and can be reloaded


## Activity 4: Dataset Integration
**Deliverable:** Integrated dataset saved in project working format

### Action 1: Prepare integration inputs
- Sub-action 1.1: Select integration variables
  - Task 1.1.1: Identify the variable(s) used to define batches/datasets for integration
    - Deliverable: Documented integration variable(s)
    - Completion criteria: Variable(s) recorded and present in metadata

### Action 2: Run integration
- Sub-action 2.1: Perform integration procedure
  - Task 2.1.1: Execute integration and store integrated representation in the dataset
    - Deliverable: Dataset containing integrated representation
    - Completion criteria: Integrated representation exists and is usable for downstream clustering

### Action 3: Save integrated output
- Sub-action 3.1: Write integrated dataset to disk
  - Task 3.1.1: Save integrated dataset
    - Deliverable: Integrated output file(s)
    - Completion criteria: Output file(s) exist and can be reloaded


## Activity 5: Clustering
**Deliverable:** Clustered dataset with cluster labels stored in metadata

### Action 1: Compute clustering inputs
- Sub-action 1.1: Generate embedding / reduced representation for clustering (if required)
  - Task 1.1.1: Compute a reduced representation used for clustering
    - Deliverable: Reduced representation stored in the dataset
    - Completion criteria: Representation exists in the dataset and matches expected dimensions

### Action 2: Cluster cells
- Sub-action 2.1: Execute clustering algorithm
  - Task 2.1.1: Run clustering and store cluster labels
    - Deliverable: Dataset with cluster labels
    - Completion criteria: Cluster labels exist in metadata and contain ≥2 unique clusters

### Action 3: Save clustered output
- Sub-action 3.1: Write clustered dataset to disk
  - Task 3.1.1: Save clustered dataset
    - Deliverable: Clustered output file(s)
    - Completion criteria: Output file(s) exist and can be reloaded with cluster labels intact


## Activity 6: Differential Gene Expression Analysis
**Deliverable:** CSV file(s) of differentially expressed genes for cluster comparisons

### Action 1: Define comparisons
- Sub-action 1.1: Specify cluster comparison strategy
  - Task 1.1.1: Define which clusters are compared and how comparisons are performed
    - Deliverable: Documented comparison specification
    - Completion criteria: Comparison specification exists and can be executed programmatically

### Action 2: Compute differential expression
- Sub-action 2.1: Run DE test(s)
  - Task 2.1.1: Perform differential expression testing for the defined comparisons
    - Deliverable: DE results in memory or intermediate form
    - Completion criteria: Results include gene identifiers and effect statistics for each comparison

### Action 3: Save DE outputs
- Sub-action 3.1: Export results to CSV
  - Task 3.1.1: Write DE tables to output directory
    - Deliverable: One or more CSV files containing DE results
    - Completion criteria: CSV file(s) exist, are non-empty, and include expected columns (gene ID/name + test statistics)

