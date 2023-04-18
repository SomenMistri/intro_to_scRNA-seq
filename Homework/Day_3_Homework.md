## Day 3 Homework
Take a look at the **scRNAseq_analysis** folder:
    
    scRNAseq_analysis                         
       ├──0_Install_R_packages_and_check.Rmd
       ├──1_Human_PDAC_tissue 
       │   │                                  
       │   ├── 1_Data_loading_and_QC_filtering
       │   ├──2_Integration_and_Clustering    
       │   └──3_Data_Visualization
       │   
       └──2_Human_PDAC_PBMC**
           │                                  
           ├── 1_Data_loading_and_QC_filtering
           │     ├──PDAC_pbmc_1_filtered_feature_bc_matrix
           │     ├──PDAC_pbmc_2_filtered_feature_bc_matrix
           │     ├──PDAC_pbmc_3_filtered_feature_bc_matrix
           │     ├──PDAC_pbmc_4_filtered_feature_bc_matrix
           │     ├──1_data_loading_and_QC_filtering.Rmd
           │     └──data_filtered.rds**
           ├──2_Integration_and_Clustering**    
           └──3_Data_Visualization
                
Open the **2_Human_PDAC_PBMC > 1_Data_loading_and_QC_filtering** folder. There you will find the **data_filtered.rds** file. Copy this file to the **2_Human_PDAC_PBMC > 2_Integration_and_Clustering** folder. Then run the **2_data_integration_and_clustering.Rmd** R markdown file.

Your homework is to replicate the data integration and clustering steps we performed today on the PBMC dataset. After running all the **chunks** of this script, please save it. This will create a **2_data_integration_and_clustering.nb.html** file in the same folder. Rename this **nb.html** file with your last name and download it to your device.

Moreover, if you run all the chunks, it should create a **"Figs"** folder in the current analysis directory. Copy all the **.png** files (n=5) to your personal device.

**Finally, email the **nb.html** file along with the five **.png** files to **smistri@uvm.edu**.**

_Note: Please USE a different **"resolution"** (0.5 or 0.8 or 1.0) value for clustering than the one we used for the PDAC human dataset._

_Note: Please USE a different **"dims"** value for clustering than the one we used for the PDAC human dataset._

_Note: Please perform SCTransform normalization **without regressing out cell cycle scores**._


**This homework is worth 50 points.**  
**For this homework you will have until 11:59 PM on Tuesday, April 25th to submit.**
