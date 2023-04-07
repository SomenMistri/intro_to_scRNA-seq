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

Your homework is to replicate the data integration and clustering steps we performed today on the PBMC dataset. After running all the **chunks** of this script, please save it. This will create a **2_data_integration_and_clustering.nb.html** file in the same folder. Rename this **.html** file with your last name and email it to the **smistri@uvm.edu** address.

_Note: Please USE a different "resolution" (0.5 or 0.8 or 1.0) than the one we used for the PDAC human dataset._

**This homework is worth 50 points.**  
**For this homework you will have until 11:59PM on Sunday, April 23rd to submit.**  
