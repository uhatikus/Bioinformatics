library(data.table)
library(magrittr)
library(XML)
library(tidyverse)
library(dplyr)
library(openxlsx)
library(magrittr)
library(readxl)

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

gene_list <- unlist(as.list(read_excel("gene_list.xlsx", col_names = FALSE)))
tcga_data <- read.table("tcga_lgg_data.txt", header = TRUE, sep = "\t")
# print(tcga_data)
diff_exp_data <- read_excel("expression_data.xlsx", col_names = TRUE)
# print(diff_exp_data)
diff_exp_data <- diff_exp_data[!is.na(diff_exp_data[, 1]), ]
diff_exp_data <- diff_exp_data[!is.na(diff_exp_data[, 2]), ]
for (cutoff_value in c(0.5, 0.6, 0.7, 0.8, 0.85, 0.9)) {
  print("cutoff:")
  print(cutoff_value)
  matched_samples <- match_exp_pattern(
    gl = gene_list, # from one of the tabs in gene_list.xlsx
    d = tcga_data, # from tcga_lgg_data.txt
    de = diff_exp_data, # from the corresponding list of expression_data.xlsx
    cutoff = cutoff_value
  ) # cutoff - the minimal fraction of genes with matching pattern within selected sample
  # print(matched_samples)

  tcga_data_matched <- tcga_data[!tcga_data$bcr_patient_barcode %in% matched_samples, ]
  row.names(tcga_data_matched) <- tcga_data_matched$bcr_patient_barcode
  all_gene_median_columns <- colnames(tcga_data_matched)[grep("median", colnames(tcga_data_matched))]
  all_gene_columns <- remove_substring_from_list(all_gene_median_columns, "median")

  useful_exp_data <- diff_exp_data[diff_exp_data$hgnc_symbol %in% all_gene_columns, ]
  useful_exp_data[useful_exp_data == -1] <- 0
  # print("Not in expression_data.xlsx:")
  not_in_exp_data <- setdiff(all_gene_columns, useful_exp_data$hgnc_symbol)
  # print(not_in_exp_data)

  # Transpose the dataframe
  useful_exp_data_good_view <- as.data.frame(t(useful_exp_data), stringsAsFactors = FALSE)
  # Set the column names
  colnames(useful_exp_data_good_view) <- useful_exp_data_good_view[1, ]
  useful_exp_data_good_view <- useful_exp_data_good_view[-1, ]
  # Reset row names
  rownames(useful_exp_data_good_view) <- NULL

  useful_exp_data_for_comparison <- useful_exp_data_good_view[rep(seq_len(nrow(useful_exp_data_good_view)), each = nrow(tcga_data_matched)), ]
  new_col_names <- paste0(colnames(useful_exp_data_for_comparison), "median")
  colnames(useful_exp_data_for_comparison) <- new_col_names
  row.names(useful_exp_data_for_comparison) <- tcga_data_matched$bcr_patient_barcode
  useful_exp_data_for_comparison[] <- lapply(useful_exp_data_for_comparison, function(x) ifelse(is.character(x), as.integer(x), x))
  tcga_data_matched_for_comparison <- tcga_data_matched[colnames(useful_exp_data_for_comparison)]

  errors <- sum(abs(tcga_data_matched_for_comparison - useful_exp_data_for_comparison))
  print("Number of people: ")
  print(nrow(tcga_data_matched))
  print("Sum errors: ")
  print(errors)
  print("Errors per person: ")
  print(errors / nrow(tcga_data_matched))
  print("Errors per gene: ")
  print(errors / length(useful_exp_data$hgnc_symbol))
  print("Probability of error: ")
  print(errors / length(useful_exp_data$hgnc_symbol) / nrow(tcga_data_matched))
}

# R t-test tcga_data_matched vs tcga_data_others
# tcga_data_matched <- tcga_data[tcga_data$bcr_patient_barcode %in% matched_samples, ]
# tcga_data_others <- tcga_data[!tcga_data$bcr_patient_barcode %in% matched_samples, ]

# all_gene_median_columns <- colnames(tcga_data_matched)[grep("median", colnames(tcga_data_matched))]

# p_values <- c()

# for (col in all_gene_median_columns) {
#   # col <- gsub("median", "", col)
#   t_from_test <- t.test(tcga_data_matched[[col]], tcga_data_others[[col]],
#     alternative = c("two.sided"),
#     mu = 0, paired = FALSE, var.equal = FALSE,
#     conf.level = 0.95
#   )
#   p_values <- c(p_values, t_from_test$p.value)
# }

# pvals <- p_values
# exp_p <- c(1:length(pvals)) / length(pvals)
# observed_p_values <- -log10(pvals)
# expected_p_values <- -log10(exp_p)
# observed_p_values <- sort(observed_p_values, decreasing = TRUE)
# plot(expected_p_values, observed_p_values)
# lines(c(0:10), c(0:10), col = "red")
# print(summary(lm(observed_p_values ~ expected_p_values)))
