library(data.table)
library(openxlsx)
library(readxl)
library(ggplot2)
library(ggpubr)
library(survival)
library(survminer)


source("src/plots.R")
source("src/utils.R")

all_sheets <-c("res_MT_None.csv", "res_MT_Progress.csv", "res_Progress_None.csv")
hist_types <- c("Oligodendroglioma", "Astrocytoma", "Oligoastrocytoma")
cutoff_values <- seq(0.5, 0.8, 0.01)

tcga_data_init <- read.table("data/tcga_lgg_data.txt", header = TRUE, sep = "\t")
cgga_data_init <- read.table("data/cgga_data.txt", header = TRUE, sep = "\t")

for (sheet in all_sheets) {
  #tcga_data <- tcga_data_init
  tcga_data <- cgga_data_init
  result_tcga_init <- prepare_diff_exp_vec_and_data_vec(data=tcga_data,sheet = sheet)
  exp_vec_tcga <- result_tcga_init[[1]]
  d_exp_tcga <- result_tcga_init[[2]]
  
  result_tcga_sort = sort_patients_by_cutoff(exp_vec_tcga, d_exp_tcga)
  errors_d_exp_tcga <- result_tcga_sort[[1]]
  sorted_d_exp_tcga <- result_tcga_sort[[2]]
  
  p_values <- c()
  n_patient <- c()
  best_p_value <- 1000
  best_cutoff_value <- 0
  best_n_patient <- 0
  for (cutoff_value in cutoff_values) {
    tcga_data_for_coxph <- prepare_data(tcga_data, cutoff_value, errors_d_exp_tcga)
    n_patient <- c(n_patient, nrow(tcga_data_for_coxph))
    res.cox <-coxph(Surv(time, alive) ~ age_at_initial_pathologic_diagnosis + gender + group, data=tcga_data_for_coxph)
    
    p_value <- summary(res.cox)$coefficients["group",5]
    p_values <- c(p_values, p_value)
    
    if (!is.na(p_value)) {
      if ((best_p_value > p_value) && (length(row.names(errors_d_exp_tcga[errors_d_exp_tcga$cutoff >= cutoff_value, ]))) > 10){
        best_p_value <- p_value
        best_cutoff_value <- cutoff_value
      }
    }
  }
  
  tcga_data_for_coxph <- prepare_data(tcga_data, best_cutoff_value, errors_d_exp_tcga)
  
  combined_plot_row_1 <- ggarrange(plot_exp_vec(exp_vec_tcga), plot_d_exp(sorted_d_exp_tcga), nrow = 1, labels = c("A", "B"), widths = c(0.16, 1-0.16), heights=c(1,1))
  combined_plot_row_2 <- ggarrange(ggplot() + theme_void(), plot_arrow_with_text_boxes(xlim = c(-0.00, 1.00)), nrow=1, ncol=2, widths = c(0.16, 1-0.16), heights=c(1,1)) 
  combined_plot_row_3 <- ggarrange(plot_survival(tcga_data_for_coxph), plot_ggforest(tcga_data_for_coxph), nrow=1, ncol=2, labels = c("С", "D"))
  combined_plot_row_4 <- ggarrange(plot_cutoff_vs_p_value(cutoff_values, p_values, n_patient, max(n_patient)), nrow=1, ncol=1, labels = c("E"))
  combined_plot <- ggarrange(combined_plot_row_1, combined_plot_row_2, combined_plot_row_3, combined_plot_row_4, ncol=1, nrow = 4, heights=c(1,0.2, 1.4, 2)) +  bgcolor("white")
  
  ggsave(paste0("plots/results/", sheet, ".png"), combined_plot, width = 18, height = 26, dpi = 300)
  
  # break
}
# TODO:
# 2. Ось cutoff  0.0 -> 0.8
# 3. пунктиры к нижней части панели
# 5. ggforest подписи красивые
# 6. add hazard ratio and p-value to surv plot (C) 
#     hazard_ratio <- summary(res.cox)$coefficients["group",2]
# machine learning (classes)? RandomForest -> gene importance 
