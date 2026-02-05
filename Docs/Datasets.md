##Source data
The CellxGene census database (https://chanzuckerberg.github.io/cellxgene-census/) will be used to build a multi-tissue mast cell atlas for healthy, adult humans from 10x v2, v3, and multiome single cell sequencing assays containing a total of ~60 million cells total from 429 individual datasets. 

Data obtained from the Cell x Gene database will be downloaded in ('.h5ad') format and will contain the following metadata features at minuimum: Biological Sex, Tissue of Origin, Developmental Stage, Cell Suspension, Assay Type, and a unique dataset ID. Each ('.h5ad') file will have a gene expression matrix of raw counts that will be used for downstream quality control, processing, and analysis. 

**Test dataset** 
For testing and development a subset of a multi-tissue immune cell atlas PMID: 40804529 will be used. This subset will contain a known number of mast cells from a high quality published dataset. The test dataset includes ~ 654 samples with approximately 1.3 million immune cells profiled across different tissues derived from 22 individual donors. An aggregate gene expression matrix from this dataset will be downloaded from the CellxGene census repository containing metadata that describes tissue of origin, biological sex, and age. A cell-level subset will be performed using a simple random sample for representative mucosal and connective tissues. This subset will be utilized to implement package features and functions.  

