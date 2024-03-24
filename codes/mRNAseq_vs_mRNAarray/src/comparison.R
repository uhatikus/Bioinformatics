library(dplyr)
library(ggplot2)

mrna_seq <- read.table("gbm_tcga/data_mrna_seq_v2_rsem.txt", header = TRUE, sep = "\t") # , row.names = 1

# TODO: Data.table -> fread
mrna_microarray <- read.table("gbm_tcga/data_mrna_affymetrix_microarray.txt", header = TRUE, sep = "\t")
# mrna_microarray <- read.table("gbm_tcga/data_mrna_agilent_microarray.txt", header = TRUE, sep = "\t") # , row.names = 1

# dupicate analysis
duplicates_mrna_seq <- mrna_seq[duplicated(mrna_seq$Hugo_Symbol) | duplicated(mrna_seq$Hugo_Symbol, fromLast = TRUE), ]
duplicates_mrna_seq_ordered <- duplicates_mrna_seq[order(duplicates_mrna_seq$Hugo_Symbol), ]

duplicates_mrna_microarray <- mrna_microarray[duplicated(mrna_microarray$Hugo_Symbol) | duplicated(mrna_microarray$Hugo_Symbol, fromLast = TRUE), ]
duplicates_mrna_microarray_ordered <- duplicates_mrna_microarray[order(duplicates_mrna_microarray$Hugo_Symbol), ]

# preparing DataFrames
mrna_seq <- mrna_seq[, !names(mrna_seq) %in% "Entrez_Gene_Id"]
missing_rows_mrna_seq <- is.na(mrna_seq$Hugo_Symbol)
mrna_seq <- mrna_seq[!missing_rows_mrna_seq, ]
# mrna_seq$Hugo_Symbol = make.unique(mrna_seq$Hugo_Symbol, sep = "_DUPL#")
mrna_seq_unique <- mrna_seq %>%
  distinct(Hugo_Symbol, .keep_all = TRUE)
rownames(mrna_seq_unique) <- mrna_seq_unique$Hugo_Symbol
mrna_seq_unique <- mrna_seq_unique[, -1]

mrna_microarray <- mrna_microarray[, !names(mrna_microarray) %in% "Entrez_Gene_Id"]
missing_rows_mrna_microarray <- is.na(mrna_microarray$Hugo_Symbol)
mrna_microarray <- mrna_microarray[!missing_rows_mrna_microarray, ]
# mrna_microarray$Hugo_Symbol = make.unique(mrna_microarray$Hugo_Symbol, sep = "_DUPL#")
mrna_microarray_unique <- mrna_microarray %>%
  distinct(Hugo_Symbol, .keep_all = TRUE)
rownames(mrna_microarray_unique) <- mrna_microarray_unique$Hugo_Symbol

# Filter columns and rows that exist in both dataframes
patients_common <- intersect(colnames(mrna_seq_unique), colnames(mrna_microarray_unique))
genes_commom <- intersect(row.names(mrna_seq_unique), row.names(mrna_microarray_unique))

mrna_seq_common <- as.data.frame(t(mrna_seq_unique[genes_commom, patients_common]), stringsAsFactors = FALSE)
mrna_microarray_common <- as.data.frame(t(mrna_microarray_unique[genes_commom, patients_common]), stringsAsFactors = FALSE)

cors <- c()

for (gene_name in genes_commom[1:1000]){
  sorted_mrna_seq_common <- mrna_seq_common[order(mrna_seq_common[[gene_name]]), ]
  sorted_mrna_microarray_common <- mrna_microarray_common[order(mrna_microarray_common[[gene_name]]), ]


  rank_mrna_seq_common <- 1:length(rownames(sorted_mrna_seq_common))
  corresponding_rank_mrna_microarray_common <- match(rownames(sorted_mrna_seq_common), rownames(sorted_mrna_microarray_common))

  cors <- c(cors, cor(rank_mrna_seq_common, corresponding_rank_mrna_microarray_common))
}

data <- data.frame(values = cors)
write.table(data, file = "cors_20230907_____.csv", sep = ",", row.names = FALSE, col.names = TRUE)
# Histogram of correlation
p <- ggplot(data, aes(x = values)) +
  geom_histogram(binwidth = 0.05, fill = "lightblue", color = "gray") +
  labs(title = "Histogram ", x = "Correlation", y = "Frequency")

ggsave("Histogram_of_correlations_____.png", plot = p, width = 6, height = 4, dpi = 300)

mclust_res <- mclust::Mclust(cors)
plot(mclust_res$BIC)
# E, V - Gauss fit params
# E - specific parameters
# V - any shape


# ggplot(data.frame('values' = cors, 'class' = mclust_res$classification), aes(x = values, fill = as.character(class))) +
#   +     geom_histogram(binwidth = 0.05, alpha = 0.5) +
#   +     labs(title = "Histogram ", x = "Correlation", y = "Frequency")


# gplot(data.frame('values' = cors, 'class' = mclust_res$classification), aes(x = values, fill = as.character(class))) +
#   +     geom_density(binwidth = 0.05, alpha = 0.5) +
#   +     labs(title = "Histogram ", x = "Correlation", y = "Frequency")