# CellStrata
BIOINF 576 Project - CellStrata

## Introduction :
Mast Cells are a heterogeneous type of immune cell important for our bodies defense against parasites and are also implicated in autoimmune and allergic diseases. Mast Cell phenotype varies based on the anatomical niche. For example, Mast Cells residing in the connective tissues of the Skin, Muscle, and Heart express the protease chymase whereas mast cells residing at the mucosal surfaces of the intestine and lung do not. Furthermore, studies suggest that there are multiple mast cell phenotypes within one particular tissue hinting at functional differences between observable phenotypes. 

Due to the rarity of human tissue specimens derived from healthy donors our understanding of mast cell phenotype in homeostasis and disease is lacking. To this end, a wrapper tool that can download data from online repositories and identify a given cell type based on key marker genes would help to better understand heterogeneity in mast cell phenotype. Such a tool could also be generalized and used to stratify other cell types in an identical fashion. 

**Questions Being Addressed**: 
What are the discreet trascriptional gene expression profiles for mast cells and/or how does tissue microenvironment influence these gene expression profiles?

**Project Description**: 
Cell Strata takes a .csv list of dataset identifiers specific to GEO / Cell Census repositories querying the GEO or Cell Census databases and extracting single cell gene expression matricies in .h5ad format. The resulting .5had file is then subject to quality control, normalization, integration via removal of batch effect, and cell type annotation to identify key marker genes that are known to segerate mast cells from other cell types. Mast cells that are identified are further clustered and subjected to differential gene expression analysis to identify discreet transcriptional states and better understand how tissue of origin impacts gene expression profiles. 

**Inputs**: Cell Strata takes a .yaml file input to query the CellxGene Census databases for single cell gene expression data. Data is downloaded in the .5had format. Resulting data is subject to further data driven filtering criteria,  and added into seperate containers within the original .5had file. 

**Final Output**
The final output will be a directory of csv files representing differentially expressed genes between the discreet clusters of mast cells that are identified and tissue wise comparison of gene expression profiles. 
