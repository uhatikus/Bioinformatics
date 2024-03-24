library(survival)
library(survminer)
library(dplyr)

sigmoid_transform <- function(x) {
  mean_value <- mean(x, na.rm = TRUE)
  std_value <- sd(x, na.rm = TRUE)
  scaled_x <- (x - mean_value) / std_value
  return(1 / (1 + exp(-scaled_x)))
}

# 1.
# Log 2 
# Quantile normalization: https://artyomovlab.wustl.edu/phantasus/phantasus-tutorial.html
# 2.
# Return to pattern matching. 
# Independence from platform
# How many genes should we have to get meaningful results?
# Add pattern matching to website

mean_HFE_mrna_seq <- mean(clinical_patient_with_target_ganes_expression$HFE_mrna_seq, na.rm = TRUE)
clinical_patient_with_target_ganes_expression$HFE_mrna_seq_binary <- ifelse(clinical_patient_with_target_ganes_expression$HFE_mrna_seq > mean_HFE_mrna_seq, 1, 0) 
clinical_patient_with_target_ganes_expression$HFE_mrna_seq_sigmoid <- sigmoid_transform(clinical_patient_with_target_ganes_expression$HFE_mrna_seq)

mean_HFE_affymetrix_microarray <- mean(clinical_patient_with_target_ganes_expression$HFE_affymetrix_microarray, na.rm = TRUE)
clinical_patient_with_target_ganes_expression$HFE_affymetrix_microarray_binary <- ifelse(clinical_patient_with_target_ganes_expression$HFE_affymetrix_microarray > mean_HFE_affymetrix_microarray, 1, 0) 
clinical_patient_with_target_ganes_expression$HFE_affymetrix_microarray_sigmoid <- sigmoid_transform(clinical_patient_with_target_ganes_expression$HFE_affymetrix_microarray)

mean_HFE_agilent_microarray <- mean(clinical_patient_with_target_ganes_expression$HFE_agilent_microarray, na.rm = TRUE)
clinical_patient_with_target_ganes_expression$HFE_agilent_microarray_binary <- ifelse(clinical_patient_with_target_ganes_expression$HFE_agilent_microarray > mean_HFE_agilent_microarray, 1, 0) 
clinical_patient_with_target_ganes_expression$HFE_agilent_microarray_sigmoid <- sigmoid_transform(clinical_patient_with_target_ganes_expression$HFE_agilent_microarray)
