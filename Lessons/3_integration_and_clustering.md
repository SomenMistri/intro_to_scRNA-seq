# Data integration and clustering

## A typical single-cell RNA-seq analysis workflow
After filtering cells using QC metrics, we need to cluster the cells based on gene expression similarity. Before clustering, we might need to perform another step "integration". We will determine whether integration is necessary or not by first clustering without integration.

 <p align="center">
<img src="../img/sc_workflow_2022.jpg" width="600">
</p>

***
_**Goals of this lesson:**_ 
 
 - _To **generate cell type-specific clusters** and use known cell type marker genes to determine the identities of the clusters._
 - _To **determine whether clusters represent true cell types or cluster due to biological or technical variation**, such as clusters of cells in the S phase of the cell cycle, clusters of specific batches, or cells with high mitochondrial content._

_**Challenges:**_
 
 - _**Identifying poor quality clusters** that may be due to uninteresting biological or technical variation_
 - _**Identifying the cell types** of each cluster_
 - _Maintaining patience as this can be a highly iterative process between clustering and marker identification (sometimes even going back to the QC filtering)_

_**Recommendations:**_
 
 - _Have a good idea of your expectations for the **cell types to be present** prior to performing the clustering. Know whether you expect cell types of low complexity or higher mitochondrial content AND whether the cells are differentiating_
 - _If you have **more than one condition**, it's often helpful to perform integration to align the cells_
 - _**Regress out** number of UMIs (by default with sctransform), mitochondrial content, and cell cycle, if needed and appropriate for experiment, so not to drive clustering_
 - _Identify any junk clusters for removal or re-visit QC filtering. Possible junk clusters could include those with high **mitochondrial content** and low UMIs/genes. If comprised of a lot of cells, then may be helpful to go back to QC to filter out, then re-integrate/cluster._
 - _If **not detecting all cell types as separate clusters**, try changing the resolution or the number of PCs used for clustering_
 
***

### Copy the input file to the analysis folder
In the "1_Data_loading_and_QC_filtering" folder, you will find the "data_filtered.rds" file. Copy that file to this folder ("2_Integration_and_Clustering").

In the "2_Integration_and_Clustering" folder, you will find the **2_data_integration_and_clustering.Rmd** R markdown file. Open the R markdown file and run each chunk sequentially.



    scRNAseq_analysis                         
       ├──0_Install_R_packages_and_check.Rmd
       ├──1_Human_PDAC_tissue 
       │   │                                  
       │   ├── 1_Data_loading_and_QC_filtering
       │   │    ├──PDAC_tissue_1_filtered_feature_bc_matrix
       │   │    ├──PDAC_tissue_2_filtered_feature_bc_matrix
       │   │    ├──PDAC_tissue_3_filtered_feature_bc_matrix
       │   │    ├──PDAC_tissue_4_filtered_feature_bc_matrix
       │   │    ├──1_data_loading_and_QC_filtering.Rmd
       │   │    ├──data_filtered.rds**
       │   │    ├──data_merged.rds
       │   │    ├──data_preliminary_clustered.rds
       │   │    └──Figures
       │   │
       │   ├──2_Integration_and_Clustering**  
       │   │   └──2_data_integration_and_clustering.Rmd  
       │   │ 
       │   └──3_Data_Visualization
       │   
       └── 2_Human_PDAC_PBMC


### Load required packages

Run chunk 1 to load the required libraries:

```r
library(XVector)
library(Seurat)
library(tidyverse)
library(Matrix)
library(RCurl)
library(scales)
library(sctransform)
```
_Note: If you have not installed the packages yet, then install them first before loading._

### Load the filtered seurat object (data_filtered.rds)
Read in the filtered data by running chunk 2. The readRDS() function used in this chunk can read in R objects that were previously created.

```r
data.filtered <- readRDS ("data_filtered.rds")
```
Let's take a look at the data to make sure everything looks good:


```r
> head(data.filtered)
                      orig.ident nCount_RNA nFeature_RNA      S.Score    G2M.Score Phase  percent.MT percent.RIBO
H1_AAACGAAGTCATAGTC-1    Human-1      72725         7990 -0.056463192 -0.087881120    G1 14.34994844    14.304572
H1_AAAGAACCATTAAAGG-1    Human-1      12734         4046 -0.020854275 -0.100676194    G1  3.68305324     8.002199
H1_AAAGGATTCGGCTTGG-1    Human-1       5653         2200 -0.015267656  0.003957696   G2M  3.29028834    19.051831
H1_AAAGGGCAGTAGCAAT-1    Human-1      20289         4772 -0.079636793 -0.120778876    G1  5.86524718     9.428754
H1_AAAGGGCAGTGAATAC-1    Human-1      50408         6714 -0.054380550 -0.092793537    G1  6.21925091    18.729170
H1_AACAACCGTTGCCTAA-1    Human-1      18025         2016 -0.001272634 -0.049945061    G1  0.03883495    34.452150
H1_AACAAGAGTGTATCCA-1    Human-1      24244         6228 -0.064505132 -0.106257050    G1  0.75482594     7.597756
H1_AACAAGATCCATGAGT-1    Human-1      14398         3920 -0.050772605 -0.050942010    G1  5.36185581    15.106265
H1_AACCAACAGCGCTTCG-1    Human-1      75830         8602 -0.074791198 -0.102485391    G1  8.76170381    16.835026
H1_AACCAACCACTGGCCA-1    Human-1      19708         4857  0.522360640  0.924381183   G2M  5.97219403    10.031459
```


## Normalization
An essential first step in most mRNA expression analyses is normalization, whereby systematic variations are adjusted to **make expression counts comparable across genes and/or samples**. The counts of mapped reads for each gene are proportional to the expression of RNA ("interesting") in addition to many other factors ("uninteresting"). Normalization is the process of adjusting raw count values to account for the "uninteresting" factors. The main factors often considered during normalization are: sequencing depth and gene length.

### Methods for scRNA-seq normalization

Various methods have been developed specifically for scRNA-seq normalization. Some **simpler methods resemble what we have seen with bulk RNA-seq**; the application of **global scale factors** adjusting for a count-depth relationship that is assumed common across all genes. However, if those assumptions are not true then this basic normalization can lead to over-correction for lowly and moderately expressed genes and, in some cases, under-normalization of highly expressed genes ([Bacher R et al, 2017](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5473255/)). **More complex methods will apply correction on a per-gene basis.** In this lesson we will explore both approaches.


**Simple transformations** are those which apply the same function to each individual measurement. Common examples include a **log transform** (which is applied in the original Seurat workflow), or a square root transform (less commonly used). In the [Hafemeister and Satija, 2019 paper](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1874-1) the authors explored the issues with simple transformations. Specifically they evaluated the standard log normalization approach and found that genes with different abundances are affected differently and that **effective normalization (using the log transform) is only observed with low/medium abundance genes. Additionally, **substantial imbalances in variance were observed with the log-normalized data**. That means **we cannot treat all genes the same.**

The proposed solution was the use of **Pearson residuals for transformation**, as implemented in Seurat's `SCTransform` function. With this approach:
* Measurements are multiplied by a gene-specific weight
* Each gene is weighted based on how much evidence there is that it is non-uniformly expressed across cells
* More evidence == more of a weight; Genes that are expressed in only a small fraction of cells will be favored (useful for finding rare cell populations)
* Not just a consideration of the expression level is, but also the distribution of expression

In this lesson, we will primarily use the "SCTransform" function of seurat to normalize the data. Note that this single _SCTransform()_ command replaces _NormalizeData(), ScaleData(), and FindVariableFeatures()_ commands of the original Seurat workflow, which performs the log-normalization.

Run **chunk 3** to normalize the filtered cells (data.filtered) using the **SCTransform()** function of seurat.

```r
data_SCT <- SCTransform(data.filtered, verbose = TRUE) 
```


## Clustering cells based on top Principal components (PCs)

### Principal Component Analysis (PCA)
Principal Component Analysis (PCA) is a technique used to emphasize variation as well as similarity, and to bring out strong patterns in a dataset; it is one of the methods used for *"dimensionality reduction"*. We will briefly go over PCA in this lesson (adapted from StatQuests/Josh Starmer's YouTube video), but we strongly encourage you to explore the video [StatQuest's video](https://www.youtube.com/watch?v=_UVHneBUBW0) for a more thorough explanation/understanding. 

#### Basic explanation with a simple example

Let's say you had quantified the expression of four genes in two samples (or cells), you could plot the expression values of those genes with one sample represented on the x-axis and the other sample on the y-axis as shown below:

<p align="center">
<img src="../img/PCA_2sample_genes.png" width="600">
</p>

You could draw a line through the data in the direction representing the **most variation**, which is on the diagonal in this example. The maximum variation in the dataset is between the genes that make up the two endpoints of this line.  

We also see the genes vary somewhat above and below the line. We could draw another line through the data representing **the second most amount of variation** in the data, since this plot is in 2D (2 axes).

The genes near the ends of each line would be those with the highest variation; these genes have the **greatest influence** on the direction of the line, mathematically. For example, a small change in the value of *Gene C* would greatly change the direction of the longer line, whereas a small change in *Gene A* or *Gene D* would have little effect on it.

<p align="center">
<img src="../img/PCA_2sample_variation2.png" width="600">
</p>


We could also rotate the entire plot and view the lines representing the variation as left-to-right and up-and-down. We see most of the variation in the data is left-to-right (longer line) and the second most variation in the data is up-and-down (shorter line). You can now think of these lines as the axes that represent the variation. These axes are essentially the "Principal Components", with PC1 representing the most variation in the data and PC2 representing the second most variation in the data. 

<p align="center">
<img src="../img/PCA_2sample_rotate.png" width="300">
</p>

Now, what if we had three samples/cells, then we would have an extra direction in which we could have variation (3D). Therefore, if we have *N* samples/cells we would have *N*-directions of variation or *N* principal components (PCs)! Once these PCs have been calculated, the PC that deals with the largest variation in the dataset is designated PC1, and the next one is designated PC2 and so on. 

Now perform principal component analysis **(PCA)** on the normalized (using the SCTranform function) dataset by running **chunk 4**.

```r
data_PCA<- RunPCA(data_SCT, npcs = 40, verbose = FALSE)
```

### Identify significant PCs
To overcome the extensive technical noise in the expression of any single gene for scRNA-seq data, Seurat assigns cells to clusters based on their PCA scores derived from the expression of the integrated most variable genes, with each PC essentially representing a “metagene” that combines information across a correlated gene set. Determining how many PCs to include in the clustering step is therefore important to ensure that we are capturing the majority of the variation, or cell types, present in our dataset.

It is useful to explore the PCs prior to deciding which PCs to include for the downstream clustering analysis.

* One way of exploring the PCs is using a heatmap to visualize the most variant genes for select PCs with the genes and cells ordered by PCA scores. The idea here is to look at the PCs and determine whether the genes driving them make sense for differentiating the different cell types.

* The elbow plot is another helpful way to determine how many PCs to use for clustering so that we are capturing the majority of the variation in the data. The elbow plot visualizes the standard deviation of each PC, and we are looking for where the standard deviations begins to plateau. Essentially, where the elbow appears is usually the threshold for identifying the majority of the variation. However, this method can be quite subjective.

Run **chunk 5** to create both the **heatmap** and **elbow plot**:

```r
# Explore heatmap of PCs
DimHeatmap(data_PCA, 
           dims = 1:10, 
           cells = 500, 
           balanced = TRUE)

# Let's draw the elbow plot using the top 40 PCs
ElbowPlot(object = data_PCA, 
          ndims = 40)
```
<p align="center">
<img src="../img/heatmap_elbow_1.png" width="900">
</p>

### Cluster the cells
Clusters of cells are obtained by grouping cells based on the similarity of their gene expression profiles. Expression profile similarity is determined via distance metrics, which often take dimensionality‐reduced representations as input. Seurat assigns cells to clusters based on their PCA scores derived from the expression of the integrated most variable genes.

While PCA determines all PCs, we can only plot two at a time. In contrast, dimensionality reduction techniques such as Uniform Manifold Approximation and Projection **(UMAP)** take the information from any number of top PCs to arrange the cells in this multidimensional space. It takes those distances in multidimensional space and plot them in two dimensions working to preserve local and global structure. In this way, the distances between cells represent similarity in expression.

For clustering the cells, **Seurat** uses a graph-based clustering approach using a **K-nearest neighbor** approach, and then attempts to partition this graph into highly interconnected ‘quasi-cliques’ or ‘communities’. The first step is to construct a K-nearest neighbor (KNN) graph based on the euclidean distance in PCA space. This is done in Seurat by using the **FindNeighbors()** function. Next, Seurat iteratively groups cells together with the goal of optimizing the standard modularity function. The **FindClusters()** function of Seurat takes care of this graph-based clustering. The **"resolution"** is an important argument at this step that sets the “granularity” of the downstream clustering and will need to be optimized for every individual experiment. For datasets of 3,000 - 5,000 cells, the resolution set between 0.4-1.4 generally yields good clustering. Increased resolution values lead to a greater number of clusters, which is often required for larger datasets.

Note: The **RunUMAP()** and **FindNeighbors()** functions will require you to input number of PCA dimentions to use for the visualization. It is suggested to use the same number of PCs as input in both cases.

Run chunk 6 to perform **RunUMAP()**, **FindNeighbors()**, **FindClusters()** functions one after another:

```r
# Runs the Uniform Manifold Approximation and Projection (UMAP) dimensional reduction technique
data_clust <- RunUMAP(data_PCA, reduction = "pca", dims = 1:30)
# Determine the K-nearest neighbor graph
data_clust <- FindNeighbors(object = data_clust, 
                                dims = 1:30)
# Perform graph based clustering
data_clust <- FindClusters(object = data_clust,
                               resolution = 0.6)

# Visualize clustered cells
DimPlot(data_clust, reduction = "umap", label = TRUE) + NoLegend()

# Save the clustered plot
ggsave(path = "Figs", filename = "Clusters.png",  height=5, width=6, units='in', dpi = 300, bg = "transparent", device='png')
```

<p align="center">
<img src="../img/Clusters.png" width="600">
</p>

#### Explore metadata
Run chunk 7 to quickly check whether "seurat_clusters" column was added to the metadata or not.

```r
head(data_clust@meta.data)
```

### Explore sources of unwanted variation

The first step in the workflow is to see if our data contains any unwanted variability. The most common biological effect that is evaluated in single-cell RNA-seq data is the effect of cell cycle on the transcriptome. Another unwanted variability can be caused by batch effect which is a technical source of variation. This step of the workflow involves exploring our data to identify which covariates we would like to correct for.

Run chunk 8 to Explore the effects of cell cycle genes:

```r
# Evaluating effects of cell cycle (Phase)
DimPlot(data_clust, group.by = "Phase", label = FALSE)

# Save the plot
ggsave(path = "Figs", filename = "CellCycle_Phase.png",  height=5, width=7, units='in', dpi = 300, bg = "transparent", device='png')
```

<p align="center">
<img src="../img/CellCycle_Phase_with_ref.png" width="900">
</p>

Run chunk 9 to explore technical sources of variation such as the Batch Effect:

```r
# Set identity classes to seurat_clusters
Idents(object = data_clust) <- "seurat_clusters"

# Explore the significance of back effect on clustering
DimPlot(data_clust, split.by = "orig.ident", label = TRUE, ncol = 2)+ NoLegend()

# Save the plot
ggsave(path = "Figs", filename = "Batch_effect.png",  height=6, width=8, units='in', dpi = 300, bg = "transparent", device='png')
```
<p align="center">
<img src="../img/Batch_effect.png" width="700">
</p>

Run chunk 10 to save the SCTransformed and clustered RDS file 

```r
saveRDS(data_clust, file= "data_clust_no_integration.rds") 
```

### To integrate or not to integrate?
Generally, we always look at our clustering without integration before deciding whether we need to perform any alignment. Do not just always perform integration because you think there might be differences - explore the data.

If cells cluster by sample, condition, batch, dataset, modality, the **integration** step can greatly improve the clustering and the downstream analyses.

As the cells we are investigating in this lesson clearly cluster by batch, performing integration is a good idea.

### Data integration 

Run **chunk 11** to perform the integration in one click. 

First, split the filtered seurat object to individual human samples and perform **SCtransform** on each sample separately. At this step, we can regress out any sources of variation such as cell cycle features, percent.MT, percent.RIBO etc. Here, we will regress out variations caused by cell cycle genes (S.Score,G2M.Score) (**Note: This particular dataset does not require cell cycle regression. However, we will perform the regression anyway to practice.)**

Following SCTransform normalization, we will perform **integration** of data using the following seurat function one after another: _SelectIntegrationFeatures(), PrepSCTIntegration(), FindIntegrationAnchors(), and IntegrateData()_

```r
#Let's split the object based on orig.ident
Idents(data.filtered) <- "orig.ident"
data.list <- SplitObject(data.filtered, split.by = "ident")

# perform SCTransform normalization with cell cycle regression
data.list <- lapply(X = data.list, FUN = function(x) {
    x <- SCTransform(x, vars.to.regress = c("S.Score","G2M.Score"), verbose = FALSE)
})

# select integration features and prep step
features <- SelectIntegrationFeatures(data.list)
data.list <- PrepSCTIntegration(
  data.list,
  anchor.features = features
)

# downstream integration steps
anchors <- FindIntegrationAnchors(
  data.list,
  normalization.method = "SCT",
  anchor.features = features
)
data.integrated <- IntegrateData(anchors, normalization.method = "SCT")

#Save the integrated file
saveRDS(data.integrated, "data.integrated.rds")
```

### Cluster the cells of integrated dataset
To perform PCA and identify significant PCs in the integrated data, please run the **chunk 12**:

```r
data_PCA2<- RunPCA(data.integrated, npcs = 40, verbose = FALSE)

# Explore heatmap of PCs
DimHeatmap(data_PCA2, 
           dims = 1:10, 
           cells = 500, 
           balanced = TRUE)

# Plot the elbow plot
ElbowPlot(object = data_PCA2, 
          ndims = 40)
```

Run **chunk 13** to perform RunUMAP(), FindNeighbors(), FindClusters() functions one after another on the integrated dataset:

```{r chunk 13}
# Runs the Uniform Manifold Approximation and Projection (UMAP) dimensional reduction technique
data_clust2 <- RunUMAP(data_PCA2, reduction = "pca", dims = 1:30)
# Determine the K-nearest neighbor graph
data_clust2 <- FindNeighbors(object = data_clust2, 
                                dims = 1:30)
# Perform graph based clustering
data_clust2 <- FindClusters(object = data_clust2,
                               resolution = 0.6)

# Visualize clustered cells
DimPlot(data_clust2, reduction = "umap", label = TRUE) + NoLegend()

# Save the clustered plot
ggsave(path = "Figs", filename = "Clusters_integrated.png",  height=5, width=6, units='in', dpi = 300, bg = "transparent", device='png')
```

<p align="center">
<img src="../img/Clusters_integrated.png" width="600">
</p>

Explore technical sources of variation such as the Batch Effect on integrated data by running **chunk 14**:

```r
# Set identity classes to seurat_clusters
Idents(object = data_clust2) <- "seurat_clusters"

# Explore the significance of back effect on clustering
DimPlot(data_clust2, split.by = "orig.ident", label = TRUE, ncol = 2)+ NoLegend()

# Save the plot
ggsave(path = "Figs", filename = "Batch_effect_integrated.png",  height=6, width=8, units='in', dpi = 300, bg = "transparent", device='png')
```
<p align="center">
<img src="../img/Batch_effect_integrated.png" width="700">
</p>

### Save the integrated and clustered cells

Run **chunk 15** to save the integrated and clustered cells as a RDS file 

_(Note: this file is important for the next lesson)_

```r
saveRDS(data_clust2, file= "data_clust_integrated.rds") 
```


                                             ---------------------------------
                                             -------END of THIS LESSON -------
                                             ---------------------------------
                                             
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
