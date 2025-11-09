# TEP-Ovarian-Cancer
Code for replicating the study "Integrative Analysis of Tumor-Educated Platelets for Stage-Specific Diagnosis, Prognosis, and Therapy in Ovarian Cancer"

### 1-Data retrieval 
Open Data directory in this repository
Create a folder named "Data" on your device
Download all data in the newly created "Data" folder
Arrange downloaded data in subdirectories/subfolders as listed in this repository

### 2-Preprocessing, differential expression analysis, and ROC analysis on GSE183635
Open source directory 
Download "Data_preprocessing_differential_expression_and_ROC_analysis_for_GSE183635.R"
"Data_preprocessing_differential_expression_and_ROC_analysis_for_GSE183635.R" with Rstudio. 
Set working directory to the location of "Data"
Run all.

### 3-Preprocessing and ROC analysis on GSE158508
Open source directory 
Download "Data_preprocessing_and_ROC_analysis_for_GSE158508_for_validation.R"
Set working directory to the location of "Data" 
Run all.

### 4-Feature Selection
Download "Feature_Selection_S1.ipynb", "Feature_Selection_S2.ipynb", "Feature_Selection_S3.ipynb", and "Feature_Selection_S4.ipynb"
Run all cells for every Notebook.

### 5-Survival analysis
Open source directory 
Set working directory to "Data" folder.
Download "Survival analysis for Candidate biomarker genes (CBGs).R"
Set working directory to the location of "Data"
Run all.

### 6-ReactomeGSA pathway enrichment analysis
Open Data directory
Navigate to "Data/Reactomegsa/Input data for reactomegsa" 
Upload gene expression matrices for the four stages: (Control_S1_counts.csv, Control_S2_counts.csv, Control_S3_counts.csv, 
Control_S4_counts.csv) separately for each stage on [https://reactome.org/gsa/home] and choose camera mode.
Turn off include disease pathways.
Upload fill in metadata information for the four stages: (Control_S1.csv, Control_S2.csv, Control_S3.csv, Control_S4.csv) through the Stage column.
First comparison group is the control group, and the second comparison group is the tested stage. 
Run and wait for results. 

## Citation
If you use this work, or any of the associated code/materials, in your own research or publication, please cite our paper. 




