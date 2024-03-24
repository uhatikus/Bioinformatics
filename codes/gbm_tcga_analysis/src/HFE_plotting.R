
hfe_columns <- grep("HFE", names(clinical_patient_with_target_ganes_expression), value = TRUE)


for (hfe_col in hfe_columns) {
  surv_formula <- as.formula(paste("Surv(time, alive) ~ age + gender + ", hfe_col))
  res.cox <-coxph(surv_formula, data=clinical_patient_with_target_ganes_expression)
  stats_forest_plot <- ggforest(res.cox, data = clinical_patient_with_target_ganes_expression)
  ggsave(paste0("plots/all_patients/", "forest_plot_", hfe_col, ".png"), stats_forest_plot, width = 10, height = 5, dpi = 300)
}

clinical_patient_with_target_ganes_expression_without_NA <- na.omit(clinical_patient_with_target_ganes_expression)
clinical_patient_with_target_ganes_expression_without_NA$HFE_all_measurements_sigmoid = (clinical_patient_with_target_ganes_expression_without_NA$HFE_mrna_seq_sigmoid +
                                                            clinical_patient_with_target_ganes_expression_without_NA$HFE_agilent_microarray_sigmoid + 
                                                            clinical_patient_with_target_ganes_expression_without_NA$HFE_agilent_microarray_sigmoid) / 3

hfe_columns_without_NA <- grep("HFE", names(clinical_patient_with_target_ganes_expression_without_NA), value = TRUE)

for (hfe_col in hfe_columns_without_NA) {
  surv_formula <- as.formula(paste("Surv(time, alive) ~ age + gender + ", hfe_col))
  res.cox <-coxph(surv_formula, data=clinical_patient_with_target_ganes_expression_without_NA)
  stats_forest_plot <- ggforest(res.cox, data = clinical_patient_with_target_ganes_expression_without_NA)
  ggsave(paste0("plots/same_patients/", "forest_plot_", hfe_col, ".png"), stats_forest_plot, width = 10, height = 5, dpi = 300)
}
  
# res.cox <-coxph(Surv(time, alive) ~ age + gender + HFE_mrna_seq , data=clinical_patient_with_target_ganes_expression)
# stats_forest_plot <- ggforest(res.cox, data = clinical_patient_with_target_ganes_expression)
# ggsave(paste0("plots/", "test_stats_forest_plot_HFE_mrna_seq_rww", ".png"), stats_forest_plot, width = 10, height = 5, dpi = 300)
# 
# res.cox <-coxph(Surv(time, alive) ~ age + gender + HFE_affymetrix_microarray, data=clinical_patient_with_target_ganes_expression)
# stats_forest_plot <- ggforest(res.cox, data = clinical_patient_with_target_ganes_expression)
# ggsave(paste0("plots/", "test_stats_forest_plot_HFE_affymetrix_microarray_raw", ".png"), stats_forest_plot, width = 10, height = 5, dpi = 300)
# 
# res.cox <-coxph(Surv(time, alive) ~ age + gender + HFE_agilent_microarray, data=clinical_patient_with_target_ganes_expression)
# stats_forest_plot <- ggforest(res.cox, data = clinical_patient_with_target_ganes_expression)
# ggsave(paste0("plots/", "test_stats_forest_plot_HFE_agilent_microarray_raw", ".png"), stats_forest_plot, width = 10, height = 5, dpi = 300)



# surv_plot <- ggsurvplot(survfit(Surv(time, alive) ~ gender, data=clinical_patient),
#                         legend.title = "",
#                         #legend.labs = c("Matched", "Unmatched"),
#                         ggtheme = theme_classic(), 
#                         palette = c("#00A7E1", "#F17720")
#                         )

#ggsave(paste0("plots/", "test_surv_plot", ".png"), surv_plot$plot, width = 10, height = 5, dpi = 300)