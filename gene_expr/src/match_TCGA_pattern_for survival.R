library(data.table)
library(magrittr)
library(XML)
library(tidyverse)
library(dplyr)
library(openxlsx)
library(magrittr)
library(readxl)
library(ggplot2)

remove_substring_from_list <- function(my_list, substring) {
  modified_list <- lapply(my_list, function(item) {
    gsub(substring, "", item)
  })
  return(modified_list %>% unlist())
}

match_exp_pattern <- function(gl = gene_list, d = xml_df,
                              de = diff_exp_data[[1]], cutoff = 0.7) {
  de <- de[de$hgnc_symbol %in% gl, ]
  exp_vec <- de[["log2FoldChange"]]
  exp_vec <- lapply(exp_vec, function(x) {
    x[which(x < 0)] <- 0
    x[which(x > 0)] <- 1
    x
  })
  exp_vec <- unlist(exp_vec)

  names(exp_vec) <- paste0(de[["hgnc_symbol"]], "median")
  exp_vec <- exp_vec[names(exp_vec) %in% colnames(d)]
  d_exp <- d[, names(exp_vec)]
  matched_samples <- apply(d_exp, MARGIN = 1, function(x) {
    check <- x == exp_vec
    total <- length(check)
    equal <- length(which(check))
    if (equal / total >= cutoff) {
      return(TRUE)
    } else {
      return(FALSE)
    }
  }) %>% unlist()
  # ss <- names(matched_samples)
  # ss <- as.numeric(ss)
  d[["bcr_patient_barcode"]][which(matched_samples)]
}

analyse_sample_data <- function(patient_data_, diff_exp_data_) {
  patient_data_ <- copy(patient_data_)
  diff_exp_data_ <- copy(diff_exp_data_)
  row.names(patient_data_) <- patient_data_$bcr_patient_barcode
  all_gene_median_columns <- colnames(patient_data_)[grep("median", colnames(patient_data_))]
  all_gene_columns <- remove_substring_from_list(all_gene_median_columns, "median")

  useful_exp_data <- diff_exp_data_[diff_exp_data_$hgnc_symbol %in% all_gene_columns, ]
  useful_exp_data[useful_exp_data == -1] <- 0
  # print("Not in expression_data.xlsx:")
  # not_in_exp_data <- setdiff(all_gene_columns, useful_exp_data$hgnc_symbol)
  # print(not_in_exp_data)

  # Transpose the dataframe
  useful_exp_data_good_view <- as.data.frame(t(useful_exp_data), stringsAsFactors = FALSE)
  # Set the column names
  colnames(useful_exp_data_good_view) <- useful_exp_data_good_view[1, ]
  useful_exp_data_good_view <- useful_exp_data_good_view[-1, ]
  # Reset row names
  rownames(useful_exp_data_good_view) <- NULL

  useful_exp_data_for_comparison <- useful_exp_data_good_view[rep(seq_len(nrow(useful_exp_data_good_view)), each = nrow(patient_data_)), ]
  new_col_names <- paste0(colnames(useful_exp_data_for_comparison), "median")
  colnames(useful_exp_data_for_comparison) <- new_col_names
  row.names(useful_exp_data_for_comparison) <- patient_data_$bcr_patient_barcode
  useful_exp_data_for_comparison[] <- lapply(useful_exp_data_for_comparison, function(x) ifelse(is.character(x), as.integer(x), x))
  data_for_comparison <- patient_data_[colnames(useful_exp_data_for_comparison)]

  errors <- abs(data_for_comparison - useful_exp_data_for_comparison)
  return(errors)
}

gene_list <- unlist(as.list(read_excel("gene_list.xlsx", col_names = FALSE)))
tcga_data <- read.table("tcga_lgg_data.txt", header = TRUE, sep = "\t")
# print(tcga_data)
diff_exp_data <- read_excel("expression_data.xlsx", col_names = TRUE, sheet = "res_Progress_None.csv")
# print(diff_exp_data)
diff_exp_data <- diff_exp_data[!is.na(diff_exp_data[, 1]), ]
diff_exp_data <- diff_exp_data[!is.na(diff_exp_data[, 2]), ]

cutoff_values <- c(0.5, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85) # , 0.9, 1

for (cutoff_value_init in cutoff_values) {
  cutoff_value_1 <- cutoff_value_init
  cutoff_value_2 <- 1 - cutoff_value_1
  matched_samples <- match_exp_pattern(
    gl = gene_list, # from one of the tabs in gene_list.xlsx
    d = tcga_data, # from tcga_lgg_data.txt
    de = diff_exp_data, # from the corresponding list of expression_data.xlsx
    cutoff = cutoff_value_1
  ) # cutoff - the minimal fraction of genes with matching pattern within selected sample
  tcga_data_matched <- tcga_data[tcga_data$bcr_patient_barcode %in% matched_samples, ]
  file_path <- paste("../coxph/given_files/lgg_matched_cutoff_res_Progress_None", cutoff_value_1, ".txt")
  list_as_text <- unlist(tcga_data_matched$bcr_patient_barcode)
  writeLines(list_as_text, file_path)

  unmatched_samples <- match_exp_pattern(
    gl = gene_list, # from one of the tabs in gene_list.xlsx
    d = tcga_data, # from tcga_lgg_data.txt
    de = diff_exp_data, # from the corresponding list of expression_data.xlsx
    cutoff = cutoff_value_2
  ) # cutoff - the minimal fraction of genes with matching pattern within selected sample
  tcga_data_unmatched <- tcga_data[!tcga_data$bcr_patient_barcode %in% unmatched_samples, ]
  file_path <- paste("../coxph/given_files/lgg_unmatched_cutoff_res_Progress_None", cutoff_value_1, ".txt")
  list_as_text <- unlist(tcga_data_unmatched$bcr_patient_barcode)
  writeLines(list_as_text, file_path)
}

# TODO:
# 1. Remove random genes and see the result.
# 2. Interferon gamma expression signature in tumors (and low grade glioma в частности).
# чем больше Т клеток -- тем больше Interferon gamma expression
