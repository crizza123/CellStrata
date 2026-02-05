# Work Breakdown Structure (WBS)

## Activity 1: CellxGene Census API query for data
### Action 1: Define query inputs
- Task 1.1: Identify dataset scope and identifiers to query
    - Deliverable: Documented list of dataset identifiers or selection 
- Task 1.2: Define required metadata/fields to retrieve

### Action 2: Retrieve data
- Task 2.1: Develop function to  
  - Taske 2.2: Identify metadata fields that are used in the repository

### Action 3: Generate code to query Cell x Gene database

**Deliverables - Activity 1** 
1. Downloaded `.h5ad` files with attached metadata for analysis
2. Visualizations for the metadata 
3. 
## Activity 2: Pre-processing of Dataset 

**Deliverables:** QC-filtered dataset stored in the project’s working format (e.g., `.h5ad` with filtered cells/genes)

### Action 1: Compute QC metrics
- SubTask 1 - Use scanpy to open data matrix
- SubTask 2 - Compute QC metrics
  - 
  -
- SubTask 3 -

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

**Deliverable** Filtered `.h5ad` files with 

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


