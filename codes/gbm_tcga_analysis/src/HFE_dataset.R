target_ganes <- c("HFE")

mrna_seq_target_genes <- as.data.frame(t(mrna_seq[rownames(mrna_seq) %in% target_ganes, ]), check.names=FALSE)
affymetrix_microarray_target_genes <- as.data.frame(t(affymetrix_microarray[rownames(affymetrix_microarray) %in% target_ganes, ]))
agilent_microarray_target_genes <- as.data.frame(t(agilent_microarray[rownames(agilent_microarray) %in% target_ganes, ]))


mrna_seq_target_genes$HFE_mrna_seq = mrna_seq_target_genes$HFE
affymetrix_microarray_target_genes$HFE_affymetrix_microarray = affymetrix_microarray_target_genes$HFE
agilent_microarray_target_genes$HFE_agilent_microarray = agilent_microarray_target_genes$HFE



clinical_patient_with_target_ganes_expression <- clinical_patient %>% mutate(Symbol = rownames(clinical_patient))
clinical_patient_with_target_ganes_expression <- full_join(clinical_patient_with_target_ganes_expression, mrna_seq_target_genes %>% 
                                                             mutate(Symbol = rownames(mrna_seq_target_genes)), by = "Symbol") %>% 
  full_join(., affymetrix_microarray_target_genes %>% 
              mutate(Symbol = rownames(affymetrix_microarray_target_genes)), by = "Symbol") %>% 
  full_join(., agilent_microarray_target_genes%>% 
              mutate(Symbol = rownames(agilent_microarray_target_genes)), by = "Symbol")
clinical_patient_with_target_ganes_expression <- clinical_patient_with_target_ganes_expression[, !names(clinical_patient_with_target_ganes_expression) %in% c("HFE", "HFE.x", "HFE.y")]
