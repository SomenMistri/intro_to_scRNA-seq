# Data visualization and interpretation

## A typical single cell RNA-seq analysis workflow
Now that we have identified our desired clusters, we can move on to marker identification, which will allow us to verify the identity of certain clusters and help surmise the identity of any unknown clusters.
 <p align="center">
<img src="../img/sc_workflow_2022.jpg" width="600">
</p>

***

_**Goals:**_ 
 
 - _To **determine the gene markers** for each of the clusters_
 - _To **identify cell types** of each cluster using markers_
 - _To determine whether there's a need to **re-cluster based on cell type markers**, perhaps clusters need to be merged or split_

_**Challenges:**_
 
 - _Over-interpretation of the results_
 - _Combining different types of marker identification_

_**Recommendations:**_
 
 - _Think of the results as hypotheses that need verification. Inflated p-values can lead to over-interpretation of results (essentially each cell is used as a replicate). Top markers are most trustworthy._
 - _Identify all markers conserved between conditions for each cluster_
 - _Identify markers that are differentially expressed between specific clusters_

***

### Copy the input file to the analysis folder
In the "2_Integration_and_Clustering" folder, you will find the "data_clust_integrated.rds" file. Copy that file to this folder ("3_Data_visualization").

In the "3_Data_visualization" folder, you will find the **3_Data_visualization.Rmd** R markdown file. Open the R markdown file and run each chunks sequentially.



    scRNAseq_analysis                         
       ├──0_Install_R_packages_and_check.Rmd
       ├──1_Human_PDAC_tissue 
       │   │                                  
       │   ├── 1_Data_loading_and_QC_filtering
       │   │ 
       │   │
       │   ├──2_Integration_and_Clustering**  
       │   │    ├──2_data_integration_and_clustering.Rmd
       │   │    ├──data_clust_integrated.rds**
       │   │    ├──data_clust_no_integration.rds
       │   │    ├──data_filtered.rds
       │   │    ├──data.integrated.rds
       │   │    ├──data_filtered.rds
       │   │    └──Figs
       │   │ 
       │   └──3_Data_Visualization**
       │        └──3_Data_Visualization.rmd
       │   
       └── 2_Human_PDAC_PBMC


### Load required packages

Run chunk 1 to load the required libraries

```r
library(XVector)
library(Seurat)
library(tidyverse)
library(Matrix)
library(RCurl)
library(scales)
library(sctransform)
```

### Load the integrated clustered seurat object (data_clust_integrated.rds) 
Read in the filtered data by running chunk 2. The readRDS() function used in this chunk can read in R objects that were previously created.

```r
data_clust_integrated <- readRDS ("data_clust_integrated.rds")

# Visualize the clustered (integrated) cells
DimPlot(data_clust_integrated, reduction = "umap", label = TRUE) + NoLegend()

# Save the plot
ggsave(path = "Figs", filename = "Clusters_integrated.png",  height=5, width=6, units='in', dpi = 300, bg = "transparent", device='png')
```
<p align="center">
<img src="../img/Clusters_integrated.png" width="600">
</p>

### Find markers for each cluster
Run chunk 3 to find markers for every cluster compared to all remaining cells. Report only the positive ones"

```r
DefaultAssay(data_clust_integrated)  <- "RNA"
Idents(data_clust_integrated) <- "seurat_clusters"
data_markers <- FindAllMarkers(data_clust_integrated, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
data_markers %>%
    group_by(cluster) %>%
    slice_max(n = 5, order_by = p_val)
```







****
### Citation

To cite material from this course in your publications, please use:

> Jihe Liu, William J. Gammerdinger, Meeta Mistry, Mary E. Piper, & Radhika S. Khetani. (2022, January 6). hbctraining/Intro-to-shell-flipped: Shell and HPC Lessons from HCBC (first release). Zenodo. https://doi.org/10.5281/zenodo.5826091

A lot of time and effort went into the preparation of these materials. Citations help us understand the needs of the community, gain recognition for our work, and attract further funding to support our teaching activities. Thank you for citing this material if it helped you in your data analysis.

---
These materials have been developed by members of the teaching team at the [Harvard Chan Bioinformatics Core (HBC)](http://bioinformatics.sph.harvard.edu/). These are open access materials distributed under the terms of the [Creative Commons Attribution license](https://creativecommons.org/licenses/by/4.0/) (CC BY 4.0), which permits unrestricted use, distribution, and reproduction in any medium, provided the original author and source are credited.

Some materials used in these lessons were derived from work that is Copyright © Data Carpentry (http://datacarpentry.org/). 
All Data Carpentry instructional material is made available under the [Creative Commons Attribution license](https://creativecommons.org/licenses/by/4.0/) (CC BY 4.0).
****