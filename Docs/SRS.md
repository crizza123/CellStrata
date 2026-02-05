**Project Description**: 

Cell Strata takes a .csv list of dataset identifiers specific to GEO / Cell Census repositories querying  Cell Census databases to obtain associated metadata in .csv and single cell gene expression matricies primarily in .h5ad format. Options for metadata filtering are also available.

The resulting .5had file is then subject to quality control, normalization, integration, and cell type annotation to identify key marker genes that are known to segerate mast cells from other cell types. Mast cells that are identified are further clustered and subjected to differential gene expression analysis to identify discreet transcriptional states and better understand how tissue of origin impacts gene expression profiles. 

**Goal of the Project**
The goal of the project is to extract, process, and identify mast cells from publicly available databases for use in better understanding discreet transcriptional states that these cells may occupy in different tissue niches.

**Features**
1. CellxGene Census API query for Data 
2. Quality Control 
4. Dataset Integration
5. Annotation across datasets
