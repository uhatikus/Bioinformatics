#install.packages("data")
#library(data.table)


mrna_seq <- fread("data/data_mrna_seq_v2_rsem.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE, data.table = FALSE)
mrna_seq <- mrna_seq[, !names(mrna_seq) %in% c("Entrez_Gene_Id")]
mrna_seq <- mrna_seq[!duplicated(mrna_seq$Hugo_Symbol) & !duplicated(mrna_seq$Hugo_Symbol, fromLast = TRUE), ]
mrna_seq <- mrna_seq[complete.cases(mrna_seq), ]
rownames(mrna_seq) = mrna_seq$Hugo_Symbol
mrna_seq <- mrna_seq[, !names(mrna_seq) %in% c("Hugo_Symbol")]
new_column_names <- gsub("\\.01$", "", names(mrna_seq))  # Remove .01 from column names
new_column_names <- gsub("\\-", ".", new_column_names)  # Replace - with .
names(mrna_seq) <- new_column_names
new_colnames <- substr(colnames(mrna_seq), 1, nchar(colnames(mrna_seq)) - 3)
colnames(mrna_seq) <- new_colnames

affymetrix_microarray <- fread("data/data_mrna_affymetrix_microarray.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE, data.table = FALSE)
affymetrix_microarray <- affymetrix_microarray[!duplicated(affymetrix_microarray$Hugo_Symbol) & !duplicated(affymetrix_microarray$Hugo_Symbol, fromLast = TRUE), ]
rownames(affymetrix_microarray) = affymetrix_microarray$Hugo_Symbol
affymetrix_microarray <- affymetrix_microarray[-c(1), !names(affymetrix_microarray) %in% c("Entrez_Gene_Id", "Hugo_Symbol")]
new_column_names <- gsub("\\.01$", "", names(affymetrix_microarray))  # Remove .01 from column names
new_column_names <- gsub("\\-", ".", new_column_names)  # Replace - with .
names(affymetrix_microarray) <- new_column_names
new_colnames <- substr(colnames(affymetrix_microarray), 1, nchar(colnames(affymetrix_microarray)) - 3)
colnames(affymetrix_microarray) <- new_colnames

agilent_microarray <- fread("data/data_mrna_agilent_microarray.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE, data.table = FALSE)
agilent_microarray <- agilent_microarray[!duplicated(agilent_microarray$Hugo_Symbol) & !duplicated(agilent_microarray$Hugo_Symbol, fromLast = TRUE), ]
rownames(agilent_microarray) = agilent_microarray$Hugo_Symbol
agilent_microarray <- agilent_microarray[-c(1), !names(agilent_microarray) %in% c("Entrez_Gene_Id", "Hugo_Symbol")]
new_column_names <- gsub("\\.01$", "", names(agilent_microarray))  # Remove .01 from column names
new_column_names <- gsub("\\.02$", "", names(agilent_microarray))  # Remove .01 from column names
new_column_names <- gsub("\\-", ".", new_column_names)  # Replace - with .
names(agilent_microarray) <- new_column_names
new_colnames <- substr(colnames(agilent_microarray), 1, nchar(colnames(agilent_microarray)) - 3)
colnames(agilent_microarray) <- new_colnames

clinical_patient <- fread("data/data_clinical_patient.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE, data.table = FALSE)
rownames(clinical_patient) = clinical_patient$PATIENT_ID
clinical_patient <- clinical_patient[, !names(clinical_patient) %in% c("OTHER_PATIENT_ID", "PATIENT_ID", "FORM_COMPLETION_DATE", "HISTORY_LGG_DX_OF_BRAIN_TISSUE", "PROSPECTIVE_COLLECTION", "RETROSPECTIVE_COLLECTION", "RACE", "ETHNICITY", "HISTORY_OTHER_MALIGNANCY", "HISTORY_NEOADJUVANT_TRTYN", "INITIAL_PATHOLOGIC_DX_YEAR", "METHOD_OF_INITIAL_SAMPLE_PROCUREMENT", "METHOD_OF_INITIAL_SAMPLE_PROCUREMENT_OTHER", "TUMOR_STATUS", "KARNOFSKY_PERFORMANCE_SCORE", "ECOG_SCORE", "PERFORMANCE_STATUS_TIMING", "RADIATION_TREATMENT_ADJUVANT", "PHARMACEUTICAL_TX_ADJUVANT", "TREATMENT_OUTCOME_FIRST_COURSE", "NEW_TUMOR_EVENT_AFTER_INITIAL_TREATMENT", "PRIMARY_SITE_PATIENT", "DISEASE_CODE", "HISTOLOGICAL_DIAGNOSIS", "ICD_10", "ICD_O_3_HISTOLOGY", "ICD_O_3_SITE", "INFORMED_CONSENT_VERIFIED", "PROJECT_CODE", "TISSUE_SOURCE_SITE", "SITE_OF_TUMOR_TISSUE", "DFS_STATUS", "DFS_MONTHS", "DAYS_TO_INITIAL_PATHOLOGIC_DIAGNOSIS")]
new_column_names <- gsub("\\-", ".", rownames(clinical_patient))  # Replace - with .
rownames(clinical_patient) <- new_column_names
clinical_patient <- clinical_patient[clinical_patient$OS_STATUS != "[Not Available]", ]
clinical_patient$time <- as.numeric(clinical_patient$OS_MONTHS)
clinical_patient$alive <- ifelse(clinical_patient$OS_STATUS == "1:DECEASED", 1, 0)
clinical_patient$gender <- ifelse(clinical_patient$SEX == "Male", 1, 2)
clinical_patient$age <- clinical_patient$AGE
clinical_patient <- clinical_patient[, !names(clinical_patient) %in% c("OS_STATUS", "OS_MONTHS", "SEX", "AGE")]
