## Day 1 Homework

### Launch a RStudio Server session on UVM VACC Open OnDemand
Log into UVM VACC Open OnDemand [website](https://vacc-ondemand.uvm.edu/). Click on the "Rstudio Server" option from the "Interactive Apps" dropdown menu. Click on launch with the default options. Your job request should be on queue for a brief moment before starting the session. Click on "Connect to Rstudio Server" icon" to start.
   
### Install required R packages and check whether the installation is a success

Install the following R packages using the the `BiocManager::install()` function:

```r
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("XVector",force = TRUE)
BiocManager::install("multtest")
BiocManager::install("glmGamPoi")
```

Now install the packages listed below from **CRAN** using the `install.packages()` function:

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

Finally, please check that all the packages were installed successfully by **loading** them using the `library()` function:

```r
library(XVector)
library(Seurat)
library(tidyverse)
library(Matrix)
library(RCurl)
library(scales)
library(sctransform)
```
If there are no errors in loading the packages, then installation was a success and we are ready for day 2.
                
This homework is entirely **optional**. This homework aims to ensure that, on **day 2**, our scripts will run smoothly. If you see error messages while installing and loading the packages, please email me **(smistri@uvm.edu)** with the error message as soon as possible. Thank you.**