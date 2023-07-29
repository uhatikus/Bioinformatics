prepare_diff_exp_vec_and_data_vec <- function(data, sheet) {
  gene_list <- unlist(as.list(read_excel("data/gene_list.xlsx", col_names = FALSE, sheet = sheet)))
  diff_exp_data <- read_excel("data/expression_data.xlsx", col_names = TRUE, sheet = sheet)
  
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
  
  exp_vec <- exp_vec[names(exp_vec) %in% colnames(data)]
  d_exp <- data[, c("bcr_patient_barcode", names(exp_vec))]
  
  rownames(d_exp) <- d_exp$bcr_patient_barcode
  d_exp$bcr_patient_barcode <- NULL
  colnames(d_exp) <- gsub("median", "", colnames(d_exp))
  
  exp_vec <- data.frame(t(exp_vec))
  colnames(exp_vec) <- gsub("median", "", colnames(exp_vec))
  
  return(list(exp_vec, d_exp))
  
}


sort_patients_by_cutoff <- function(exp_vec, d_exp){
  errors_d_exp <- abs(d_exp - exp_vec[rep(seq_len(nrow(exp_vec)), each = nrow(d_exp)), ])
  errors_d_exp$cutoff <- 1 - rowSums(errors_d_exp)/ncol(errors_d_exp)
  right_order <- order(errors_d_exp$cutoff, decreasing = TRUE)
  errors_d_exp <- errors_d_exp[right_order, ]
  sorted_d_exp <- d_exp[right_order, ]
  return(list(errors_d_exp, sorted_d_exp))
}

prepare_data <- function(data, cutoff_value, errors_d_exp){
  data_for_coxph <- data
  rownames(data_for_coxph) <- data_for_coxph$bcr_patient_barcode
  data_for_coxph$group <- 0  # Initialize flag column with 0 or sample(c(-1, 0, 1), nrow(data_for_coxph), replace = TRUE)
  data_for_coxph$group[row.names(data_for_coxph) %in% row.names(errors_d_exp[errors_d_exp$cutoff >= cutoff_value, ])] <- -1  # Set -1 if row name is in list1
  data_for_coxph$group[row.names(data_for_coxph) %in% row.names(errors_d_exp[errors_d_exp$cutoff <= 1-cutoff_value, ])] <- 1  # Set 1 if row name is in list2
  data_for_coxph$gender <- ifelse(data_for_coxph$gender == "MALE", 1, 2)
  data_for_coxph <- data_for_coxph[data_for_coxph$group != 0, ]
  return(data_for_coxph)
}