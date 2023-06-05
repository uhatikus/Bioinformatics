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
diff_exp_data <- read_excel("expression_data.xlsx", col_names = TRUE)
# print(diff_exp_data)
diff_exp_data <- diff_exp_data[!is.na(diff_exp_data[, 1]), ]
diff_exp_data <- diff_exp_data[!is.na(diff_exp_data[, 2]), ]

cutoff_values <- c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.65, 0.7, 0.71, 0.72, 0.73, 0.74, 0.75, 0.8, 0.85, 0.9, 1)
p_errors_matched <- c()
p_errors_unmatched <- c()
cutoff_values <- c(0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.65, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85) # , 0.9, 1

# cutoff_values <- c(0.7, 0.8)
n_genes_unmatched <- c()

for (cutoff_value in cutoff_values) {
  # print("cutoff:")
  # print(cutoff_value)
  matched_samples <- match_exp_pattern(
    gl = gene_list, # from one of the tabs in gene_list.xlsx
    d = tcga_data, # from tcga_lgg_data.txt
    de = diff_exp_data, # from the corresponding list of expression_data.xlsx
    cutoff = cutoff_value
  ) # cutoff - the minimal fraction of genes with matching pattern within selected sample
  tcga_data_matched <- tcga_data[tcga_data$bcr_patient_barcode %in% matched_samples, ]
  errors_matched <- analyse_sample_data(tcga_data_matched, diff_exp_data)
  p_error <- sum(errors_matched) / 80 / nrow(tcga_data_matched)
  p_errors_matched <- c(p_errors_matched, p_error)

  print(cutoff_value)
  strange_genes <- colnames(errors_matched[colSums(errors_matched) > 0.5 * nrow(errors_matched)])
  if (length(strange_genes) < 20) {
    print(length(strange_genes))
    print(gsub("median", "", strange_genes))
  }
  strange_people <- rownames(errors_matched[rowSums(errors_matched) > 0.5 * ncol(errors_matched), ])
  print(length(strange_people) == 0)

  tcga_data_unmatched <- tcga_data[!tcga_data$bcr_patient_barcode %in% matched_samples, ]
  errors_unmatched <- analyse_sample_data(tcga_data_unmatched, diff_exp_data)
  p_error <- sum(errors_unmatched) / 80 / nrow(tcga_data_unmatched)
  p_errors_unmatched <- c(p_errors_unmatched, p_error)
  n_genes_unmatched <- c(n_genes_unmatched, length(strange_genes))
}
# p <- ggplot(data = data.frame(cutoff_values, p_errors_matched), aes(x = cutoff_values, y = p_errors_matched)) +
#   geom_point(color = "blue") +
#   geom_smooth(method = "lm", se = FALSE, color = "red", linetype = "dashed") +
#   labs(title = "Scatter Plot with Trend Line", x = "X", y = "Y") +
#   theme_minimal()
p <- ggplot() +
  geom_point(data = data.frame(x = cutoff_values, y = n_genes_unmatched), aes(x = cutoff_values, y = n_genes_unmatched), color = "blue") +
  # geom_point(data = data.frame(x = cutoff_values, y = p_errors_matched), aes(x = cutoff_values, y = p_errors_matched), color = "blue") +
  # geom_point(data = data.frame(x = cutoff_values, y = p_errors_unmatched), aes(x = cutoff_values, y = p_errors_unmatched), color = "red") +
  # geom_smooth(data = data.frame(x = cutoff_values, y = p_errors_matched), aes(x = cutoff_values, y = p_errors_matched), method = "lm", se = FALSE, color = "red", linetype = "dashed") +
  # geom_smooth(data = data.frame(x = cutoff_values, y = p_errors_unmatched), aes(x = cutoff_values, y = p_errors_unmatched), method = "lm", se = FALSE, color = "purple", linetype = "dashed") +
  # labs(title = "Probability of Error vs Cut-off", x = "Cut-off value", y = "Probability of Error") +
  labs(title = "Number of not matched genes vs Cut-off", x = "Cut-off value", y = "Number of not macthed genes ") +
  # coord_cartesian(ylim = c(0, 1)) +
  theme_minimal()
ggsave("plot.png", p, width = 8, height = 6, dpi = 300)
# IRdisplay::display(p)

# plotting: ggsci, ggpubr, ggplot -- recommended

# 1. survival analysis cutoff 0.8 vs cutoff 0.2, 0.7 vs 0.3, 0.6 vs 0.4, etc.
# 2. remove genes from expression_data and repeat 1.
# slides -> googlelides, 00:00 from Mon to Tue, English
# - expain
# - patter exp
# - tcga data
# - matched vs anti-matched
# - cutoff vs number of genes that couldn't match
# - same genes? - yes -> list of genes (gene -- number of unmatched)
# - survival analysis 1 and 2.

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
