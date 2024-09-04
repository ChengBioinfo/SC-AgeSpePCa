# SC-AgeSpePCa
Single-cell analysis identify the age-associated divergent microenvironments and progression patterns in prostate cancer

The scripts involved in this study are available as followed (each R file has a detailed bibliography that can be identified by Rstudio):

1.data_processing.R

Creating Seurat objects from raw matrix, and carrying out quality control.

2.harmony.R

Removing batch effects using Harmony. Clustering, annotation and functional analysis of all cells.

3.Doublet and InferCNV.R

Removing doublets and infering malignant epithelia with inferCNV.

4.1.Epi_NMF.R

Deciphering the tumor heterogeneity with NMF. Batch effect removing, clustering, annotation, trajectory analysis, functional analysis, and all other analyses concerning malignant epithelia.

4.2.Myeloid_harmony.R

Batch effect removing, clustering, annotation, trajectory analysis, functional analysis, and all other analyses concerning myeloid cells.

4.3.Tcell_harmony.R

Batch effect removing, clustering, annotation, trajectory analysis, functional analysis, and all other analyses concerning T and NK cells.

4.4.Fibro_harmony.R

Batch effect removing, clustering, annotation, trajectory analysis, functional analysis, and all other analyses concerning mesenchymal cells.

5.1.CellCommunication.R

Cell communication analyses between tumor cells and macrophages, macrophages and T cells, etc.

5.2.CellCommunication_epi_all.R

Cell communication analyses inferring how microenvironment cells regulating tumor cells, through NicheNet.

6.ST.R

Clustering, annotation, and all other analysis concerning spatial transcriptomes.
