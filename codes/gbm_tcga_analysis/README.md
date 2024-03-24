# Code for Analysis of the HFE gene in GBM (Glioblastoma)

### How to use

#### Prepare the data

You should have the following files in the `data` folder (downloaded from [CBioPortal](https://www.cbioportal.org/study/summary?id=gbm_tcga_pub2013)):

- gbm_tcga_analysis/data/data_clinical_patient.txt
- gbm_tcga_analysis/data/data_mrna_affymetrix_microarray.txt
- gbm_tcga_analysis/data/data_mrna_agilent_microarray.txt
- gbm_tcga_analysis/data/data_mrna_seq_v2_rsem.txt

#### Run the code

You should run the R code in the following order:

1. `src/prepare_data.R`
2. `src/HFE_dataset.R`
3. `src/HFE_analysis.R`
4. `src/HFE_plotting.R`

Plots will appear in folder `plots`
