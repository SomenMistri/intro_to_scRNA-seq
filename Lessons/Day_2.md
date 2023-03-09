# Day 2

## A typical single cell RNA-seq analysis workflow
After quantifying gene expression from raw sequencing reads, we need to bring the count matrix data (filtered_feature_bc_matrix) into R to generate metrics for performing QC and further downstream analysis.

 <p align="center">
<img src="../img/sc_workflow_2022.jpg" width="600">
</p>


## Filter cells using quality metrics

### Launching a RStudio Server session on UVM VACC Open OnDemand
Log into UVM VACC Open OnDemand [website](https://vacc-ondemand.uvm.edu/). Click on the "Rstudio Server" option from the "Interactive Apps" dropdown menu. Click on launch with the default options. Your job request should be on queue for a brief moment before starting the session. Click on "Connect to Rstudio Server" icon" to start.

On your lower right hand side, you should see the files of your home directory. Please click the "scRNAseq_analysis" folder icon. You should see the following folder structure:

    scRNAseq_analysis                         
       ├──0_Install_R_packages_and_check.Rmd
       ├──1_Human_PDAC_tissue                                   
       │   ├── 1_QC_filtering_and_clustering
       │   │    ├──PDAC_tissue_1_filtered_feature_bc_matrix
       │   │    ├──PDAC_tissue_2_filtered_feature_bc_matrix
       │   │    ├──PDAC_tissue_3_filtered_feature_bc_matrix
       │   │    ├──PDAC_tissue_4_filtered_feature_bc_matrix
       │   │    └──1_QC_filtering_and_first_pass_clustering.Rmd
       │   ├──2_integration_and_clustering    
       │   └──3_data_visualization
       │   
       └── 2_Human_PDAC_PBMC
     
Please navigate to the "/scRNAseq_analysis/1_Human_PDAC_tissue/" folder and click on the 1_QC_filtering_and_first_pass_clustering.Rmd file. The R markdown file can also be found [here](https://github.com/SomenMistri/intro_to_scRNA-seq/blob/main/scripts/1_QC_filtering_and_first_pass_clustering.Rmd).
    
               
### Load data and make seurat objects

To load the required packages using the **library()** function, Run **chunk 1** by clicking on the "Run Current Chunk" button on the right. This will load the following packages. 
```
library(XVector)
library(Seurat)
library(tidyverse)
library(Matrix)
library(RCurl)
library(scales)
library(sctransform)
```
_Note: If you have not installed the packages yet, then install them first before loading_
     
First, let's read in the data by running **chunk 2**. The Read10X() function reads in the output of the cellranger pipeline from 10X, returning a unique molecular identified (UMI) count matrix. The values in this matrix represent the number of molecules for each feature (i.e. gene; row) that are detected in each cell (column).
```{r chunk 2}
data1 <- Read10X("PDAC_tissue_1_filtered_feature_bc_matrix")
data2 <- Read10X("PDAC_tissue_2_filtered_feature_bc_matrix")
data3 <- Read10X("PDAC_tissue_3_filtered_feature_bc_matrix")
data4 <- Read10X("PDAC_tissue_4_filtered_feature_bc_matrix")

#Let's take a look at how the data looks like
head(data1)
```

Then, let's use the count matrix to create a Seurat object by running **chunk 3**. The seurat object serves as a container for both the data (like the count matrix) and analysis (e.g. PCA, metadata) for a single-cell dataset.
```
data_seurat1 <- CreateSeuratObject(counts = data1, project = "Human-1", min.cells = 3, min.features = 200)
data_seurat2 <- CreateSeuratObject(counts = data2, project = "Human-2", min.cells = 3, min.features = 200)
data_seurat3 <- CreateSeuratObject(counts = data3, project = "Human-3", min.cells = 3, min.features = 200)
data_seurat4 <- CreateSeuratObject(counts = data4, project = "Human-4", min.cells = 3, min.features = 200)
```

### (Optional) Cell Cycle Scoring
In some cases, there is a need for mitigating the effects of cell cycle heterogeneity in scRNA-seq data.This can be done by calculating cell cycle phase scores based on known cell cycle markers , and regressing these out of the data during pre-processing.
 <p align="center">
<img src="../img/cellcycle.png" width="500">
</p>

To perform cell cycle scoring, run **chunk 4**. In this chunk, we are first Log Normalizing individual seurat objects using the **NormalizeData()** function. Then, we are using the **CellCycleScoring()** function to assign each cell a cell cycle score, based on its expression of G2/M and S phase markers. Seurat stores the s.genes and g2m.genes in the "cc.genes.updated.2019" list. Finally, we are using the **head()** function to make sure that **Phase** information was added as a new column. 
```
#segregate the "cc.genes.updated.2019" list into markers of G2/M phase and markers of S phase
s.genes <- cc.genes.updated.2019$s.genes
g2m.genes <- cc.genes.updated.2019$g2m.genes

#Prior to running "CellCycleScoring" command, each seurat object needs to be Lognormalized using "NormalizeData" function
data_norm1 <- NormalizeData(data_seurat1, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)
data_norm2 <- NormalizeData(data_seurat2, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)
data_norm3 <- NormalizeData(data_seurat3, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)
data_norm4 <- NormalizeData(data_seurat4, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)

#Now perform CellCycleScoring for each seurat objects
data_norm1 <- CellCycleScoring(data_norm1, s.features = s.genes, g2m.features = g2m.genes, set.ident = FALSE, verbose = FALSE)
data_norm2 <- CellCycleScoring(data_norm2, s.features = s.genes, g2m.features = g2m.genes, set.ident = FALSE, verbose = FALSE)
data_norm3 <- CellCycleScoring(data_norm3, s.features = s.genes, g2m.features = g2m.genes, set.ident = FALSE, verbose = FALSE)
data_norm4 <- CellCycleScoring(data_norm4, s.features = s.genes, g2m.features = g2m.genes, set.ident = FALSE, verbose = FALSE)

# view cell cycle scores and phase assignments
head(data_norm1[[]])
```

### Merge seurat objects and generate quality metrics

Run **chunk 5** to merge all four seurat objects into one. The merge() function merges the raw count matrices of two or more Seurat objects creating a new Seurat object with a combined raw count matrix. Then, let's take a look at the metadata of the merged seurat object using the View() function.


```
#NOTE: By default, merge() function combines Seurat objects based on the raw count matrices, erasing any previous normalization
data_merged <- merge(data_norm1, y = c(data_norm2, data_norm3, data_norm4), add.cell.ids = c("H1", "H2", "H3","H4"), project = "Human_1234")

# Explore merged metadata
View(data_merged@meta.data)
```

