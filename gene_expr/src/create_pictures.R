library(data.table)
library(magrittr)
library(XML)
library(tidyverse)
library(dplyr)
library(openxlsx)
library(magrittr)
library(readxl)
library(ggplot2)
library(reshape2)
library(ggpubr)
library(survival)
library(survminer)


all_sheets <-c("res_MT_None.csv", "res_MT_Progress.csv", "res_Progress_None.csv")
hist_types <- c("Oligodendroglioma", "Astrocytoma", "Oligoastrocytoma")
# Grade II oligodendrogliomas are low grade tumors.
# astrocytomas range from grade 1 (most benign) to grade 4 (most malignant). (both?)
# oligoastrocytomas as grade II (low-grade) or grade III (anaplastic)
cutoff_values <- seq(0.5, 0.8, 0.01)

tcga_data_init <- read.table("tcga_lgg_data.txt", header = TRUE, sep = "\t")


prepare_diff_exp_vec_and_tcga_data_vec <- function(tcga_data, sheet) {
  gene_list <- unlist(as.list(read_excel("gene_list.xlsx", col_names = FALSE, sheet = sheet)))
  diff_exp_data <- read_excel("expression_data.xlsx", col_names = TRUE, sheet = sheet)
  
  diff_exp_data <- diff_exp_data[!is.na(diff_exp_data[, 1]), ]
  diff_exp_data <- diff_exp_data[!is.na(diff_exp_data[, 2]), ]
  
  diff_exp_data <- diff_exp_data[diff_exp_data$hgnc_symbol %in% gene_list, ]
  
  exp_vec <- diff_exp_data[["log2FoldChange"]]
  exp_vec <- lapply(exp_vec, function(x) {
    x[which(x < 0)] <- 0
    x[which(x > 0)] <- 1
    x
  })
  exp_vec <- unlist(exp_vec)
  
  names(exp_vec) <- paste0(diff_exp_data[["hgnc_symbol"]], "median")
  
  exp_vec <- exp_vec[names(exp_vec) %in% colnames(tcga_data)]
  d_exp <- tcga_data[, c("bcr_patient_barcode", names(exp_vec))]
  
  rownames(d_exp) <- d_exp$bcr_patient_barcode
  d_exp$bcr_patient_barcode <- NULL
  colnames(d_exp) <- gsub("median", "", colnames(d_exp))
  
  #d_exp <- head(d_exp)
  
  
  exp_vec <- data.frame(t(exp_vec))
  colnames(exp_vec) <- gsub("median", "", colnames(exp_vec))
  
  return(list(exp_vec, d_exp))
  
}


plot_exp_vec <- function(exp_vec, sheet){
  g <- ggplot(melt(as.matrix(-2*exp_vec+1)), aes(x = "", y = Var2, fill = as.factor(value))) +
    geom_tile(color = "white",
              lwd = 0.2,
              linetype = 1) +
    labs(title = "Pattern",
         x = "",
         y = "Genes") + 
    scale_fill_manual(values = c("-1" = "green", "1" = "red"), labels=c("-1" = "LogFC > 0", "1" = "LogFC < 0"), name="") +
    coord_fixed(1) +
    #theme_minimal() +
    theme(axis.text.x=element_blank(), legend.position="left", plot.title= element_text(hjust = 0.5)) 
  #ggsave(paste0(sheet, "_pattern.png"), g, width = 5, height = 5, dpi = 300)
  return(g)
}

plot_d_exp <- function(d_exp, sheet){
  g <- ggplot(melt(as.matrix(- 2*d_exp + 1)), aes(x = Var1, y = Var2, fill = as.factor(value))) +
    geom_tile(color = "white",
              lwd = 0.05,
              linetype = 1) +
    labs(title = "Patients",
         x = "",
         y = "Genes")  +
    scale_fill_manual(values = c("-1" = "#00C957", "1" = "#DC143C"), labels=c("-1" = "Gene Expression > Median", "1" = "Gene Expression < Median"), name="") +
    
    # scale_fill_gradient(low ="red", high = "green") + 
    coord_fixed(5) + 
    #theme_minimal() +
    theme(axis.text.x = element_blank(), legend.position="right", plot.title= element_text(hjust = 0.5), axis.text.y=element_blank(), axis.title.y=element_blank())
  #ggsave(paste0(sheet, "_patients.png"), g, width = 15, height = 5, dpi = 300)
  return(g)
}

plot_erros_d_exp <- function(errors_d_exp, sheet){
  errors_d_exp <- errors_d_exp[, c(ncol(errors_d_exp), 1:(ncol(errors_d_exp)-1))]
  g <- ggplot(melt(as.matrix(errors_d_exp)), aes(x = Var1, y = Var2, fill = value)) +
    geom_tile(color = "white",
              lwd = 0.1,
              linetype = 1) +
    labs(title = "Errors",
         x = "Patients",
         y = "Genes") + scale_fill_gradient(low ="#00A7E1", high = "#F17720") + 
    coord_fixed(1)+
    theme(axis.text.x = element_text(angle = 90)) + 
    geom_text(aes(label = round(value, 2)), color = "white", size = 1)
  #ggsave(paste0(sheet, "_erros.png"), g, width = 49, height = 5, dpi = 300)
  return(g)
}

plot_arrow <- function(xlim = c(-0.05, 1.05)){
  arrow_df <- data.frame(
    x = c(0.0, 1.0),
    y = c(0, 0)
  )
  
  # Create the ggplot object
  g <- ggplot(arrow_df, aes(x, y))
  
  # Add the arrow
  g <- g + geom_path(arrow = arrow(length = unit(0.3, "cm"))) + coord_cartesian(xlim = xlim, ylim = c(-0.1, 0.25))
  
  # add rectangle with text
  g <- g + geom_rect(aes(xmin = 0.0, xmax = 0.2, ymin = 0.05, ymax = 0.15),
                     fill = "#00A7E1", color = "black", alpha = 0.7)
  g <- g + geom_text(aes(x = 0.1, y = 0.1, label = "Matched"))
  
  
  # add rectangle with text
  g <- g + geom_rect(aes(xmin = 0.8, xmax = 1.0, ymin = 0.05, ymax = 0.15),
                     fill = "#F17720", color = "black", alpha = 0.7)
  g <- g + geom_text(aes(x = 0.9, y = 0.1, label = "Unmatched"))
  
  
  g <- g + geom_text(aes(x = 0.05, y = 0.25, label = "0.2"))
  g <- g + geom_text(aes(x = 0.5, y = 0.25, label = "0.5"))
  
  
  g <- g + theme_void()
  return(g)
}

sort_patients_by_cutoff <- function(exp_vec, d_exp){
  errors_d_exp <- abs(d_exp - exp_vec[rep(seq_len(nrow(exp_vec)), each = nrow(d_exp)), ])
  errors_d_exp$cutoff <- 1 - rowSums(errors_d_exp)/ncol(errors_d_exp)
  right_order <- order(errors_d_exp$cutoff, decreasing = TRUE)
  errors_d_exp <- errors_d_exp[right_order, ]
  sorted_d_exp <- d_exp[right_order, ]
  return(list(errors_d_exp, sorted_d_exp))
  
}

prepare_data <- function(tcga_data, cutoff_value){
  tcga_data_for_coxph <- tcga_data
  rownames(tcga_data_for_coxph) <- tcga_data_for_coxph$bcr_patient_barcode
  tcga_data_for_coxph$group <- 0  # Initialize flag column with 0
  tcga_data_for_coxph$group[row.names(tcga_data_for_coxph) %in% row.names(errors_d_exp[errors_d_exp$cutoff >= cutoff_value, ])] <- -1  # Set -1 if row name is in list1
  tcga_data_for_coxph$group[row.names(tcga_data_for_coxph) %in% row.names(errors_d_exp[errors_d_exp$cutoff <= 1-cutoff_value, ])] <- 1  # Set 1 if row name is in list2
  tcga_data_for_coxph$gender <- ifelse(tcga_data_for_coxph$gender == "MALE", 1, 2)
  tcga_data_for_coxph$alive <- ifelse(tcga_data_for_coxph$vital_status == "Dead", 1, 0)
  #tcga_data_for_coxph$age_at_initial_pathologic_diagnosis
  #tcga_data_for_coxph$time <- rowSums(tcga_data_for_coxph[, c("days_to_death", "days_to_last_followup")], na.rm = TRUE)
  #tcga_data_for_coxph$
  tcga_data_for_coxph <- tcga_data_for_coxph[tcga_data_for_coxph$group != 0, ]
  return(tcga_data_for_coxph)
}

#perform_coxph_analysis <- function(errors_d_exp){
  #row.names(errors_d_exp[errors_d_exp$cutoff > 0.85, ])
#}
for (sheet in all_sheets) {
  for (hist_type in hist_types) {
    tcga_data <- tcga_data_init[tcga_data_init$histological_type == hist_type, ]
    result_init <- prepare_diff_exp_vec_and_tcga_data_vec(tcga_data=tcga_data,sheet = sheet)
    exp_vec <- result_init[[1]]
    d_exp <- result_init[[2]]
    
    result_sort = sort_patients_by_cutoff(exp_vec, d_exp)
    errors_d_exp <- result_sort[[1]]
    sorted_d_exp <- result_sort[[2]]
    
    #plot_exp_vec(exp_vec, sheet)  
    #plot_d_exp(sorted_d_exp, sheet)
    
  
    
    p_values <- c()
    n_patient <- c()
    best_p_value <- 1000
    best_cutoff_value <- 0
    best_n_patient <- 0
    for (cutoff_value in cutoff_values) {
      tcga_data_for_coxph <- prepare_data(tcga_data, cutoff_value)
      n_patient <- c(n_patient, nrow(tcga_data_for_coxph))
      res.cox <-coxph(Surv(time, alive) ~ age_at_initial_pathologic_diagnosis + gender + group, data=tcga_data_for_coxph)
      #res.cox <-coxph(Surv(time, alive) ~ gender + group, data=tcga_data_for_coxph)
      #print(summary(res.cox))
      p_value <- (summary(res.cox)$coefficients[3,5])
      hazard_ratio <- summary(res.cox)$coefficients[3,2]
      p_values <- c(p_values, p_value)
      
      if (!is.na(p_value)) {
        if ((best_p_value > p_value) && (length(row.names(errors_d_exp[errors_d_exp$cutoff >= cutoff_value, ]))) > 10){
            best_p_value <- p_value
            best_cutoff_value <- cutoff_value
        }
      }
    }
    
    #best_cutoff_value = 0.7
    
    tcga_data_for_coxph <- prepare_data(tcga_data, best_cutoff_value)
    
    surv_plot <- ggsurvplot(survfit(Surv(time, alive) ~ group, data=tcga_data_for_coxph),
                     legend.title = "",
                     legend.labs = c("Matched", "Unmatched"),
                     #title = paste0("Survival analysis with cutoff value = ", best_cutoff_value, ", p_value = ", signif(best_p_value, 3)),
                     ggtheme = theme_classic(), palette = c("#00A7E1", "#F17720"))
    #print(surv_plot)
    coeff <- max(n_patient)
    cutoff_vs_p_value <- ggplot(data = data.frame(x = 1-cutoff_values, y1 = p_values, y2 = n_patient/coeff), aes(x = x)) +
      geom_point(aes(y = y1, color = "Y1")) +
      geom_point(aes(y = y2, color = "Y2")) +
      geom_line(linetype = "dashed", aes(y = y1), color = "#00A7E1") +
      geom_line(linetype = "dashed", aes(y = y2), color = "#F17720") +
      labs(title = "Two Sets of Y Data",
             x = "X-axis Label",
             y = "Y-axis Label") +
      scale_color_manual(name = "",
                           values = c("Y1" = "#00A7E1", "Y2" = "#F17720"),
                           labels = c("p-value", "Number of patients")) + 
      #geom_line(linetype = "dashed", data = data.frame(x = cutoff_values, y = p_values), aes(x = cutoff_values, y = p_values), color="#00A7E1") +
      labs(title = "", x = "level of dissimilarity", y = "COXPH p-values") + theme_bw() + theme_classic() +
      geom_hline(yintercept = 0.05, linetype = "dashed", color = "red")  + theme(plot.margin = unit(c(0,3, 1, 1), "cm")) +
      geom_text(aes(x = 0.22, y = 0.08, label = "p-value = 0.05")) + 
      #geom_point(data = data.frame(x = cutoff_values, y = n_patient/coeff), aes(x = x, y = y), color = "#F17720") + 
      #geom_line(linetype = "dashed", data = data.frame(x = cutoff_values, y = n_patient/coeff), aes(x = x, y = y), color = "#F17720")+
      scale_y_continuous(
        
        # Features of the first axis
        name = "COXPH p-values",
        
        # Add a second axis and specify its features
        sec.axis = sec_axis(~.*coeff, name="Number of patients in analysis")
      ) +
      theme(legend.position = c(0.7, 0.9))
    
    cox_model <-coxph(Surv(time, alive) ~ age_at_initial_pathologic_diagnosis + gender + group, data=tcga_data_for_coxph)
    stats_forest_plot <- ggforest(cox_model, data = tcga_data_for_coxph)
    
    combined_plot_row_1 <- ggarrange(plot_exp_vec(exp_vec, sheet), plot_d_exp(sorted_d_exp, sheet), nrow = 1, labels = c("A", "B"), widths = c(0.16, 1-0.16), heights=c(1,1))
    ifelse((sheet == "res_MT_None.csv"), xlim_arrow <- (c(-0.05, 1.05)),
        ifelse((sheet == "res_MT_Progress.csv"), xlim_arrow <- (c(-0.00, 1.00)),
          ifelse((sheet == "res_Progress_None.csv"), xlim_arrow <- (c(-0.00, 1.00)),
                 xlim_arrow <- c(-0.00, 1.00))))
    combined_plot_row_2 <- ggarrange(plot_arrow, plot_arrow(xlim = xlim_arrow), nrow=1, ncol=2, widths = c(0.16, 1-0.16), heights=c(1,1)) 
    combined_plot_row_3 <- ggarrange(surv_plot$plot+ theme(plot.margin = unit(c(0,1, 1, 2), "cm"))+ theme(legend.position = c(0.9, 0.9)), cutoff_vs_p_value, nrow=1, ncol=2, labels = c("С", "D"))
    plot_row_4 <- stats_forest_plot
    combined_plot <- ggarrange(combined_plot_row_1, combined_plot_row_2, combined_plot_row_3, plot_row_4, ncol=1, nrow = 4, heights=c(1,0.2, 1.4, 2)) +  bgcolor("white")
    
    #print(combined_plot)
    if (sheet == "res_MT_Progress.csv"){
      ggsave(paste0(sheet, "_", hist_type, ".png"), combined_plot, width = 25, height = 26, dpi = 300)
    }else{
      ggsave(paste0(sheet, "_", hist_type, ".png"), combined_plot, width = 18, height = 26, dpi = 300)
    }
    #plot_erros_d_exp(errors_d_exp, sheet)
    
  }
}



# TODO:
# 1. Combine 2 datasets (plot just CGG - to check) (median - for eachh dataset different)
# 2. Ось cutoff  0.0 -> 0.8
# 3. пунктиры к нижней части панели
# 4. E <-> D
# 5. ggforest подписи красивые


# TODO:
# 1. put legend inside the plot
# 2. legend for plot D



# TODO:
# 1. line matched unmatched
#arrow_df <- data.frame(x = c(0, 1),y = c(0, 0))
# Create the ggplot object
#p <- ggplot(arrow_df, aes(x, y))

# Add the arrow
#p + geom_path(arrow = arrow(length = unit(0.3, "cm"))) + theme_void()

# 2. p-value = 0.05
# 3. sizes 
# 4. colors
# 5. number of patient by cutoff in graph D
# 6. cutof -> scale from 1 to -1