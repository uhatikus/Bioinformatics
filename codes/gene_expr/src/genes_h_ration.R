library(data.table)
library(openxlsx)
library(readxl)
library(ggplot2)
library(survival)

all_sheets <-c("res_MT_None.csv", "res_MT_Progress.csv", "res_Progress_None.csv")

tcga_data_init <- read.table("data/tcga_lgg_data.txt", header = TRUE, sep = "\t")
cgga_data_init <- read.table("data/cgga_data.txt", header = TRUE, sep = "\t")

gene_list <- c()
for (sheet in all_sheets) {
  gene_list <- c(gene_list, unlist(as.list(read_excel("data/gene_list.xlsx", col_names = FALSE, sheet = sheet))))
}
gene_list <- unique(gene_list)

hr_tcga <- c()
not_in_tcga <- c()
in_tcga <- c()
for (g in gene_list) {
  gm <- paste0(g, "median")
  if (!gm %in% colnames(tcga_data_init)){
    not_in_tcga <- c(not_in_tcga, gm)
  } else {
    in_tcga <- c(in_tcga, gm)
    f1 <- as.formula(paste0("Surv(time, alive) ~ gender + ", gm))
    # f1 <- as.formula(paste0("Surv(time, alive) ~ age_at_initial_pathologic_diagnosis + gender + ", gm))
    cox_model <-coxph(f1, data=tcga_data_init)
    hr <- summary(cox_model)$coefficients[gm, 2]
    hr_tcga <- c(hr_tcga, hr)
  }
}

tcga_genes_hrs <- data.frame(gene = in_tcga, hazard_ratio_tcga = hr_tcga)
write.csv(tcga_genes_hrs, file = "tcga_genes_hrs_2.csv", row.names = FALSE)



hr_cgga <- c()
not_in_cgga <- c()
in_cgga <- c()
for (g in gene_list) {
  gm <- paste0(g, "median")
  if (!gm %in% colnames(cgga_data_init)){
    not_in_cgga <- c(not_in_cgga, gm)
  } else {
    in_cgga <- c(in_cgga, gm)
    f1 <- as.formula(paste0("Surv(time, alive) ~ age_at_initial_pathologic_diagnosis + gender + ", gm))
    cox_model <-coxph(f1, data=cgga_data_init)
    hr <- summary(cox_model)$coefficients[gm, 2]
    hr_cgga <- c(hr_cgga, hr)
  }
}

cgga_genes_hrs <- data.frame(gene = in_cgga, hazard_ratio_cgga = hr_cgga)
write.csv(cgga_genes_hrs, file = "cgga_genes_hrs_2.csv", row.names = FALSE)

combined_data <- merge(tcga_genes_hrs, cgga_genes_hrs, by = "gene", all=TRUE)

p <- ggplot(combined_data, aes(x = gene)) +
  geom_bar(aes(y = hazard_ratio_tcga, color = "tcga", alpha = 0.7), stat = "identity", position = "dodge", alpha = 0.7) +
  geom_bar(aes(y = hazard_ratio_cgga, color = "ccga", alpha = 0.7), stat = "identity", position = "dodge", alpha = 0.7) +
  labs(title = "tcga vs ccga",
       x = "gene",
       y = "hazard ratio") +
  scale_color_manual(name = "Datasets",
                     values = c("tcga" = "blue", "ccga" = "red")) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
print(p)

ggsave(paste0("plots/results/tcga_vs_ccga_2.png"), p, width = 10, height = 7, dpi = 300)