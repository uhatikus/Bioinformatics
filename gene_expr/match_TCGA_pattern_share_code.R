library(data.table)
library(magrittr)
library(XML)
library(tidyverse)
library(dplyr)
library(openxlsx)
library(magrittr)
library(readxl)

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
matched_samples <- match_exp_pattern(
  gl = gene_list, # from one of the tabs in gene_list.xlsx
  d = tcga_data, # from tcga_lgg_data.txt
  de = diff_exp_data, # from the corresponding list of expression_data.xlsx
  cutoff = 0.7
) # cutoff - the minimal fraction of genes with matching pattern within selected sample
print(matched_samples)

all_gene_median_columns <- colnames(tcga_data_matched)[grep("median", colnames(tcga_data_matched))]

tcga_data_matched <- tcga_data[tcga_data$bcr_patient_barcode %in% matched_samples, ]
tcga_data_others <- tcga_data[!tcga_data$bcr_patient_barcode %in% matched_samples, ]
n1 <- length(tcga_data_matched)
n2 <- length(tcga_data_others)

p_values <- c()

for (col in all_gene_median_columns) {
  # col <- gsub("median", "", col)
  t_from_test <- t.test(tcga_data_matched[[col]], tcga_data_others[[col]],
    alternative = c("two.sided"),
    mu = 0, paired = FALSE, var.equal = FALSE,
    conf.level = 0.95
  )
  p_values <- c(p_values, t_from_test$p.value)
}

pvals <- p_values
exp_p <- c(1:length(pvals)) / length(pvals)
observed_p_values <- -log10(pvals)
expected_p_values <- -log10(exp_p)
observed_p_values <- sort(observed_p_values, decreasing = TRUE)
plot(expected_p_values, observed_p_values)
lines(c(0:10), c(0:10), col = "red")
print(summary(lm(observed_p_values ~ expected_p_values)))
