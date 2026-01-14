# CellStrata
BIOINF 576 Project - CellStrata

**Introduction**:
Mast Cells are a heterogeneous type of immune cell important for host defense against parasites and also implicated in autoimmune and allergic diseases. Mast Cell phenotype varies based on the anatomical niche. For example Mast Cells residing in the connective tissues of the Skin, Muscle, and Heart express the protease chymase whereas mast cells residing at the mucosal surfaces of the intestine and lung do not. Furthermore, studies suggest that multiple   

**Question Being Addressed**:   How does a given cell type vary in terms of gene expression based on its tissue of origin? 

**Project Description**: Cell Strata takes a .csv list of dataset identifiers specific to GEO / Cell Census repositories as an input from multiple tissue types to query the GEO or Cell Census databases to extract single cell gene expression matricies in .h5ad format. The resulting .5had file is then subject to quality control, normalization, integration, and cell type annotation to identify key marker genes for a cell type of interest. Cells that are identified are clustered and subjected to differential gene expression analysis to identify how tissue of origin impacts gene expression profiles. 

**Inputs**: Cell Strata takes a .csv file input to query the GEO / Cell Census databases for single cell gene expression data downloaded in the .5had format. Resulting data is the processed and added into a seperate layer within the original .5had file. 

**Final Output**
