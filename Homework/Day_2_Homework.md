## Day 2 Homework
First, take a look at the **scRNAseq_analysis** folder:
    
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
       └──2_Human_PDAC_PBMC**
           │                                  
           └── 1_Data_loading_and_QC_filtering
                ├──PDAC_pbmc_1_filtered_feature_bc_matrix
                ├──PDAC_pbmc_2_filtered_feature_bc_matrix
                ├──PDAC_pbmc_3_filtered_feature_bc_matrix
                ├──PDAC_pbmc_4_filtered_feature_bc_matrix
                └──1_data_loading_and_QC_filtering.Rmd
                
Open the **2_Human_PDAC_PBMC** folder. There you will find four count matrices of PBMC samples from four human PDAC patients.

Also, you will find the **1_data_loading_and_QC_filtering.Rmd** R markdown file. This is the same R markdown file we went over in this lesson.

Your homework is to replicate the data loading and QC filtering steps we performed today on these PBMC dataset. After running all the **chunks** of this script, please save it. This will create a **1_data_loading_and_QC_filtering.nb.html** file in the same folder. Rename this **.html** file with your last name and email it to the **smistri@uvm.edu** address.

_Hint: You will need to edit **chunk 2** to load the new PBMC count matrices which have different filenames._

_Note: Please **DO NOT** use the same filtering parameters we used for the PDAC human dataset. Use a different set of filtering parameters for the human PBMC datasets_

**This homework is worth 50 points.**  
**For this homework you will have until 11:59PM on Sunday, April 16th to submit.**  
