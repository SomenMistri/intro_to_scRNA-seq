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
Open the 1_Install_R_packages_and_check.Rmd R markdown file. Run the first chunk (chunk 1) to install R packages Bioconductor using the the BiocManager::install() function.


     


  