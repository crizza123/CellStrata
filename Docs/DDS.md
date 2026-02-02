config / dataset list (.csv)
        |
        v
[1] data_acquisition
    - query CellxGene Census
    - obtain metadata fields and ontology terms  (.csv)
    - materialize as AnnData (.h5ad)
        |
        v
[2] qc
    - compute QC metrics
    - filter cells/genes/doublets
        |
        v
[3] normalization
    - normalize / log-transform (store in layers or X)
        |
        v
[4] integration
    - integrate across datasets/batches (store integrated embedding)
        |
        v
[5] clustering
    - build neighbors/graph on integrated space
    - assign cluster labels (obs["cluster"])
        |
        v
[6] differential_expression
    - run DE between clusters (and optionally tissue-wise comparisons)
        |
        v
[7] outputs
    - write CSVs of DE results (+ any summary tables)
