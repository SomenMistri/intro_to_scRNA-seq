# Day 2

## A typical single cell RNA-seq analysis workflow
After quantifying gene expression from raw sequencing reads, we need to bring the count matrix data (filtered_feature_bc_matrix) into R to generate metrics for performing QC and further downstream analysis.

 <p align="center">
<img src="../img/sc_workflow_2022.jpg" width="600">
</p>

### A typical "Cell Ranger count" output looks like this:
    PDAC_human_10x                                      
        └── outs                                          
            ├── filtered_feature_bc_matrix     #this is the folder you need for downstream R analysis
            │   ├── barcodes.tsv.gz
            │   ├── features.tsv.gz
            │   └── matrix.mtx.gz
            ├── filtered_feature_bc_matrix.h5  #Alternatively, this file can also be used as a input for downstream R analysis
            ├── raw_feature_bc_matrix        
            └── raw_feature_bc_matrix.h5

## Filter cells using quality metrics
            
### Getting ready with the count matrix files

### Launching a RStudio Server session on UVM VACC Open OnDemand

### Install required R packages and check whether the installation is a success

Open the 1_Install_R_packages_and_check.Rmd R markdown file.

> **Note 1:  All the package names listed below are case sensitive!**
 
> **Note 2**: At any point (especially if you’ve used R/Bioconductor in the past), in the console **R may ask you if you want to update any old packages by asking Update all/some/none? [a/s/n]:**. If you see this, **type "a" at the prompt and hit Enter** to update any old packages. _Updating packages can sometimes take quite a bit of time to run, so please account for that before you start with these installations._  

> **Note 3:** If you see a message in your console along the lines of “binary version available but the source version is later”, followed by a question, **“Do you want to install from sources the package which needs compilation? y/n”, type n for no, and hit enter**.

 Run the first chunk (**chunk 1**) to install R packages Bioconductor using the the `BiocManager::install()` function.

`if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("XVector",force = TRUE)
BiocManager::install("multtest")
BiocManager::install("glmGamPo")`

Now run the second chunk (**chunk 2**) to install packages listed below from **CRAN** using the `install.packages()` function. 
`install.packages('tidyverse')
install.packages('Matrix')
install.packages('RCurl')
install.packages('scales')
install.packages('metap')
install.packages('Seurat')
install.packages("ggplot2")`

Finally, please check that all the packages were installed successfully by **loading** them using the `library()` function (**chunk 3**).
`library(XVector)
library(Seurat)
library(tidyverse)
library(Matrix)
library(RCurl)
library(scales)`






     


  