library(ggplot2)

# gbm_data <- RTCGAToolbox::getFirehoseData(dataset = "LGG",
#                                           rundate = "20160128",
#                                           clinical = TRUE, 
#                                           Mutation = TRUE,
#                                           RNASeq2GeneNorm = TRUE, # normalized number of transcriped (log_2 transformed)
#                                           RNASeq2Gene = TRUE, # use
#                                           mRNAArray=TRUE)     # use 

rnaseq <- RTCGAToolbox::getData(gbm_data, "RNASeq2Gene")
rnaarray <- RTCGAToolbox::getData(gbm_data, "mRNAArray", platform=1)
rnaaray_exp2 <- 2^rnaarray
rnaaray_exp2 <- rnaaray_exp2[order(rownames(rnaaray_exp2)), ]

rnaseq <- as.data.frame(rnaseq, stringsAsFactors = FALSE)
rnaaray_exp2 <- as.data.frame(rnaaray_exp2, stringsAsFactors = FALSE)

# Filter columns and rows that exist in both dataframes
common_patients <- intersect(colnames(rnaseq), colnames(rnaaray_exp2))
common_genes <- intersect(row.names(rnaseq), row.names(rnaaray_exp2))

rnaseq_common <- as.data.frame(t(rnaseq[common_genes, common_patients]), stringsAsFactors = FALSE)
rnaaray_exp2_common <- as.data.frame(t(rnaaray_exp2[common_genes, common_patients]), stringsAsFactors = FALSE)



# Sort the gene row by expression value

cors <- c()

# for (gene_name in common_genes){
gene_name <- "ALKBH4"
  sorted_rnaseq <- rnaseq_common[order(rnaseq_common[[gene_name]]), ]
  sorted_rnaaray_exp2 <- rnaaray_exp2_common[order(rnaaray_exp2_common[[gene_name]]), ]
  
  
  rank_rnaseq <- 1:length(rownames(sorted_rnaseq))
  corresponding_rank_rnaaray_exp2 <- match(rownames(sorted_rnaseq), rownames(sorted_rnaaray_exp2))
  
  cors <- c(cors, cor(rank_rnaseq, corresponding_rank_rnaaray_exp2))
  
# }
# Plot 
plot_data <- data.frame(
  rank_X = rank_rnaseq,
  rank_Y = corresponding_rank_rnaaray_exp2
)


ggplot(plot_data, aes(x = rank_X, y = rank_Y)) +
  geom_point(color = "blue") +
  labs(x = "Rank mRNAseq", y = "Rank mRNA array", title=paste0("Gene: ", gene_name, "; correlation coef = ", cors[1])) +
  theme_minimal()
  # geom_point(data=plot_data, aes(x = rank_X+1, y = rank_Y+1), color = "red")

# same people
# by people rank
# X: rank in mRNAseq, Y: rank in mRNA_array

# + download from gliovis and compare 


# TCGA - main, later: CGGA, ICGC
# all cancer datasets : LGG, ...
# example of genes: cor = 0.99, mean, 0.0, <0 (min):                                                      4 figures
# Histogram of correlation (how many Gaussians? Clusters of genes? Gtex (level of expression???) ),       20 figures: for each cancer: 2 columns by 20 histograms
# Distance to 1.0 correlation across all cancer types: all gene in the end have low correlation           1 figure: Histogram: distance (residuals) and for all cancer types
# Relation to gene length                                                                                 1 figure: scatter plot + trend line (correlation). X: length, Y: distance
# Bayesian Meta-Analysis

# proof of concept by showing "wrong" paper

#       LGG  Cancer2  Cancer3
# Gene1
# Gene2
# Gene3
# Gene4
# 
# ML algorithm: clustering gene
# Gaussian Mixture model (autogmm)


# https://www.youtube.com/watch?v=JId304dp3tc