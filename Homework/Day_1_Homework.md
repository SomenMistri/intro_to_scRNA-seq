## Day 1 Homework (Optional)

### Getting ready with the count matrix and R markdown files
First, log in to your VACC account. Copy the scRNAseq_analysis.tar file from the shared folder to your home directory. This .tar file contains the count matrix files as well as the necessary R markdown files (.Rmd) required for this tutorial.

You can also find this .tar file [here](../Data/scRNAseq_analysis.tar/)

```bash
cp -r /gpfs1/cl/mmg232/course_materials/scRNAseq_analysis.tar
```

Then, extract the .tar file in your home directory. You should see a new folder named "scRNAseq_analysis" in your home directory. 

```bash
tar -xvf scRNAseq_analysis.tar
```

### Launching a RStudio Server session on UVM VACC Open OnDemand
Log into UVM VACC Open OnDemand [website](https://vacc-ondemand.uvm.edu/). Click on the "Rstudio Server" option from the "Interactive Apps" dropdown menu. Click on launch with the default options. Your job request should be on queue for a brief moment before starting the session. Click on "Connect to Rstudio Server" icon" to start.

On your lower right hand side, you should see the files of your home directory. Please click the "scRNAseq_analysis" folder icon. You should see the following folder structure:

    scRNAseq_analysis                         
       ├──0_Install_R_packages_and_check.Rmd
       ├──1_Human_PDAC_tissue 
       │   │                                  
       │   ├── 1_Data_loading_and_QC_filtering
       │   │    ├──PDAC_tissue_1_filtered_feature_bc_matrix
       │   │    ├──PDAC_tissue_2_filtered_feature_bc_matrix
       │   │    ├──PDAC_tissue_3_filtered_feature_bc_matrix
       │   │    ├──PDAC_tissue_4_filtered_feature_bc_matrix
       │   │    └──1_data_loading_and_QC_filtering.Rmd
       │   │
       │   ├──2_Integration_and_Clustering    
       │   └──3_Data_Visualization
       │   
       └── 2_Human_PDAC_PBMC
         
               


### Install required R packages and check whether the installation is a success

Please navigate to the "/scRNAseq_analysis/0_QC_filtering_and_clustering" folder and click on the 0_Install_R_packages_and_check.Rmd file. The R markdown file can also be found [here](https://github.com/SomenMistri/intro_to_scRNA-seq/blob/main/scripts/1_Install_R_packages_and_check.Rmd).

You can also install these packages by copying and pasting the code chunks into Rstudio console.

> **Note 1:**  All the package names listed in the R markdown file are case sensitive!**
 
> **Note 2:** At any point (especially if you’ve used R/Bioconductor in the past), in the console **R may ask you if you want to update any old packages by asking Update all/some/none? [a/s/n]:**. If you see this, **type "a" at the prompt and hit Enter** to update any old packages. _Updating packages can sometimes take quite a bit of time to run, so please account for that before you start with these installations._  

> **Note 3:** If you see a message in your console along the lines of “binary version available but the source version is later”, followed by a question, **“Do you want to install from sources the package which needs compilation? y/n”, type n for no, and hit enter**.

Run the first chunk (**chunk 1**) to install R packages Bioconductor using the the `BiocManager::install()` function.

```r
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("XVector",force = TRUE)
BiocManager::install("multtest")
BiocManager::install("glmGamPoi")
```

Now run the second chunk (**chunk 2**) to install packages listed below from **CRAN** using the `install.packages()` function. 

```r
install.packages('tidyverse')
install.packages('Matrix')
install.packages('RCurl')
install.packages('scales')
install.packages('metap')
install.packages('Seurat')
install.packages("ggplot2")
install.packages("sctransform")
```

Finally, please check that all the packages were installed successfully by **loading** them using the `library()` function (**chunk 3**).

```r
library(XVector)
library(Seurat)
library(tidyverse)
library(Matrix)
library(RCurl)
library(scales)
library(sctransform)
```
If there are no errors in loading the packages, then installation was a success and we are ready for day 2. You can now close this R markdown file and end the  RStudio session. 
                
This homework is entirely **optional**. This homework aims to ensure that, on **day 2**, our scripts will run smoothly. If you see error messages while installing and loading the packages, please email me **(smistri@uvm.edu)** with the error message as soon as possible. Thank you.**