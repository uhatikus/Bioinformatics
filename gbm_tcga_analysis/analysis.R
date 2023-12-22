library(survival)
library(survminer)
library(dplyr)

target_ganes <- c("HFE")


mrna_seq_target_genes <- as.data.frame(t(mrna_seq[rownames(mrna_seq) %in% target_ganes, ]), check.names=FALSE)
affymetrix_microarray_target_genes <- as.data.frame(t(affymetrix_microarray[rownames(affymetrix_microarray) %in% target_ganes, ]))
agilent_microarray_target_genes <- as.data.frame(t(agilent_microarray[rownames(agilent_microarray) %in% target_ganes, ]))

mrna_seq_target_genes$HFE_mrna_seq = mrna_seq_target_genes$HFE
affymetrix_microarray_target_genes$HFE_affymetrix_microarray = affymetrix_microarray_target_genes$HFE
agilent_microarray_target_genes$HFE_agilent_microarray = agilent_microarray_target_genes$HFE


#print(mrna_seq_target_genes[duplicated(rownames(mrna_seq_target_genes)) | duplicated(rownames(mrna_seq_target_genes), fromLast = TRUE), ])
#print(affymetrix_microarray_target_genes[duplicated(rownames(affymetrix_microarray_target_genes)) | duplicated(rownames(affymetrix_microarray_target_genes), fromLast = TRUE), ])
#print(agilent_microarray_target_genes[duplicated(rownames(agilent_microarray_target_genes)) | duplicated(rownames(agilent_microarray_target_genes), fromLast = TRUE), ])

#mrna_seq_target_genes_indices <- rownames(mrna_seq_target_genes)

#mrna_seq_target_genes_no_duplicates <- mrna_seq_target_genes[!duplicated(rownames(mrna_seq_target_genes)) & !duplicated(rownames(mrna_seq_target_genes), fromLast = TRUE), , drop = FALSE]
#mrna_seq_target_genes_no_duplicates$Original_Index <- mrna_seq_target_genes_indices[!duplicated(rownames(mrna_seq_target_genes)) & !duplicated(rownames(mrna_seq_target_genes), fromLast = TRUE)]
#mrna_seq_target_genes_no_duplicates <- as.data.frame(mrna_seq_target_genes_no_duplicates)
#mrna_seq_target_genes_no_duplicates <- mrna_seq_target_genes_no_duplicates[, !names(mrna_seq_target_genes_no_duplicates) %in% c("Original_Index")]


#mrna_seq_target_genes_no_duplicates <- mrna_seq_target_genes[!duplicated(rownames(mrna_seq_target_genes)) & !duplicated(rownames(mrna_seq_target_genes), fromLast = TRUE), ]

#common_indices_mrna_seq_affymetrix_microarray <- intersect(row.names(mrna_seq_target_genes), row.names(affymetrix_microarray_target_genes))
#common_indices_mrna_seq_agilent_microarray <- intersect(row.names(mrna_seq_target_genes), row.names(agilent_microarray_target_genes))
#common_indices_affymetrix_microarray_agilent_microarray <- intersect(row.names(affymetrix_microarray_target_genes), row.names(agilent_microarray_target_genes))
clinical_patient_with_target_ganes_expression <- clinical_patient %>% mutate(Symbol = rownames(clinical_patient))
clinical_patient_with_target_ganes_expression <- full_join(clinical_patient_with_target_ganes_expression, mrna_seq_target_genes %>% 
                                  mutate(Symbol = rownames(mrna_seq_target_genes)), by = "Symbol") %>% 
                                    full_join(., affymetrix_microarray_target_genes %>% 
                                                mutate(Symbol = rownames(affymetrix_microarray_target_genes)), by = "Symbol") %>% 
                                     full_join(., agilent_microarray_target_genes%>% 
                                                 mutate(Symbol = rownames(agilent_microarray_target_genes)), by = "Symbol")
clinical_patient_with_target_ganes_expression <- clinical_patient_with_target_ganes_expression[, !names(clinical_patient_with_target_ganes_expression) %in% c("HFE", "HFE.x", "HFE.y")]

#clinical_patient_with_target_ganes_expression <- full_join(clinical_patient_with_target_ganes_expression, agilent_microarray_target_genes)

#clinical_patient_with_target_ganes_expression <- clinical_patient_with_target_ganes_expression[, !names(clinical_patient_with_target_ganes_expression) %in% c("Original_Index")]


res.cox <-coxph(Surv(time, alive) ~ age + gender + HFE_mrna_seq + HFE_affymetrix_microarray + HFE_agilent_microarray, data=clinical_patient_with_target_ganes_expression)


surv_plot <- ggsurvplot(survfit(Surv(time, alive) ~ gender, data=clinical_patient),
                        legend.title = "",
                        #legend.labs = c("Matched", "Unmatched"),
                        ggtheme = theme_classic(), 
                        palette = c("#00A7E1", "#F17720")
                        )

ggsave(paste0("plots/", "test_surv_plot", ".png"), surv_plot$plot, width = 10, height = 5, dpi = 300)


res.cox <-coxph(Surv(time, alive) ~ age + gender + HFE_mrna_seq , data=clinical_patient_with_target_ganes_expression)
stats_forest_plot <- ggforest(res.cox, data = clinical_patient_with_target_ganes_expression)
ggsave(paste0("plots/", "test_stats_forest_plot_1", ".png"), stats_forest_plot, width = 10, height = 5, dpi = 300)

res.cox <-coxph(Surv(time, alive) ~ age + gender + HFE_affymetrix_microarray, data=clinical_patient_with_target_ganes_expression)
stats_forest_plot <- ggforest(res.cox, data = clinical_patient_with_target_ganes_expression)
ggsave(paste0("plots/", "test_stats_forest_plot_2", ".png"), stats_forest_plot, width = 10, height = 5, dpi = 300)

res.cox <-coxph(Surv(time, alive) ~ age + gender + HFE_agilent_microarray, data=clinical_patient_with_target_ganes_expression)
stats_forest_plot <- ggforest(res.cox, data = clinical_patient_with_target_ganes_expression)
ggsave(paste0("plots/", "test_stats_forest_plot_3", ".png"), stats_forest_plot, width = 10, height = 5, dpi = 300)
