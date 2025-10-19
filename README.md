## About the Project
This project focuses on gene expression profiling of endometriotic tissues across different subtypes.  
The research question this project aims to address is: *"Characterizing subtype-specific gene expression signatures in endometriotic tissues."*

The workflow includes data loading, data inspection, data filtering, normalization steps, followed by statistical analysis.
Multiple plots are provided to display data structure and results.



## Repository Structure
**main branch:** contains all finalized scripts (1-3), as well as the a merged script, along with the 'LICENSE' and '.gitignore' files. 

**Load_and_process_Script1:** one branch for only the first part of the project (Data loading, inspection etc.)

**Data_Filter_Script2:** one branch for the second part of the project (Data filtering, probes, class imbalance, replicates)

**Analysis_Script3:** one branch for only the analysis part

**merged_script:** contains a single file that combines all three parts into one script, enabling a more streamlined execution of the entire analysis

**script_crafting branch:** includes most of the R and Quarto scripts, as well as the intermediate commits and updates that were accumulated over the course of the project



## About the Data
The dataset originates from a published study and is based on Illumina HumanHT-12 microarray bulk expression data.  
It includes multiple metadata variables describing gene characteristics, which were excluded in the beginning of the analysis.

**The dataset**
- consists of a total number of 47,324 genes and 72 samples
- represent raw expression values with associated detection scores
- contains sample replicates and sometimes multiple probes per gene
- shows class imbalance between endometriotic subtypes
- does not include healthy control samples


**Data source:** 
[GSE141549 – Gene Expression Omnibus (GEO)](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE141549)  

**File used in this analysis:**  
'GSE141549_Non-normalized_data_GA_illumina_expression_platform_HumanHT-12.xls.gz'

**Original publication:** 
Gabriel M., Fey V., Heinosalo T., Adhikari P., Rytkönen K., Komulainen T., Huhtinen K., Laajala T. D., Siitari H., Virkki A., Suvitie P., Kujari H., Aittokallio T., Perheentupa A., & Poutanen M. (2020).  
*A relational database to identify differentially expressed genes in the endometrium and endometriosis lesions.*  
**Scientific Data**, 7, Article 294. [https://doi.org/10.1038/s41597-020-00623-x](https://doi.org/10.1038/s41597-020-00623-x)



## Requirements/Versions
- **R:** 4.5.1 (June 2025)
- **RStudio:** 2025.09.1+401
- **Quarto:** 1.7.32
- **Script format:** Quarto ('.qmd') with HTML output
    (the reposatory also contains helper R scripts, but the final analysis was converted in '.qmd' files)
- **R packages:** 'readxl', 'tidyverse', 'BiocManager', 'limma', 'stringr'



## About the Script
The workflow is divided into three main sections:

1) **Setup & Data Loading**
   - Check/install/load packages
   - Load dataset
   - Exclude matadata variables and  define matrix with raw expression data and associated detection scores

2) **Exploration & Preprocessing**
   - **Detection scores:** inspect density and rescale scores into a 0–1 “p-value range”
   - **Expression values:** quality control by using PCA, hierarchical clustering, and visualization of expression distribution across samples
   - **Filtering & normalization:** test different p-values and fraction cutoffs, select significant genes, apply log2-transformation and quantile normalization
   - **Probe & replicate handling:** aggregate probes, inspect replicate relationships, re-check expression distributions, visualization of class imbalance

3) **Main Analysis**
   - PCA on normalized data using all subtypes
   - Reclassification of subtypes into PE vs NonPE
   - PCA on normalized data using new classes
   - Differential Expressed Genes analysis (DEG)
   - PCA on only significant DEGs



**How to use the script:** Set the working directory to the correct data path, then run the script completely to generate all plots and results.
