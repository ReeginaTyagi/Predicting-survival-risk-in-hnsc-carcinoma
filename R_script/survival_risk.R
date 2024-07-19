# Load necessary libraries
library(TCGAbiolinks)
library(biomaRt)
library(dplyr)
library(SummarizedExperiment)
library(DESeq2)
library(maftools)
library(survival)
library(survminer)
library(glmnet)
library(Hmisc)

# Set working directory and load data
setwd("C:/Users/tyagi/Downloads/Predicting-survival-risk-in-hnsc-carcinoma")
df <- readRDS("C:/Users/tyagi/Downloads/Predicting-survival-risk-in-hnsc-carcinoma/tcga_hnsc_data.RDS")

# Prepare data
colData <- colData(df)
rowData <- rowData(df)

# Filter out normal tissue samples and prepare survival data
colData <- colData[colData$sample_type != "Solid Tissue Normal", ]
colData$survival_time <- coalesce(colData$days_to_death, colData$days_to_last_follow_up)
colData <- colData[!is.na(colData$survival_time), ]
colData$vital_status <- ifelse(colData$vital_status == "Alive", 0, 1)
colData <- colData[, c("vital_status", "survival_time")]
colData <- as.data.frame(colData)

# Annotate genes and filter for lncRNAs
genes <- rowData$external_gene_name
mart <- useMart("ENSEMBL_MART_ENSEMBL", host = "useast.ensembl.org")
mart <- useDataset("hsapiens_gene_ensembl", mart)
annotLookup <- getBM(mart = mart, attributes = c("hgnc_symbol", "entrezgene_id", "ensembl_gene_id", "gene_biotype"), filter = "hgnc_symbol", values = genes, uniqueRows = TRUE)
lncRNA <- annotLookup[annotLookup$gene_biotype == "lncRNA", ]
ensembl <- lncRNA$ensembl_gene_id
rowData <- rowData[rowData$ensembl_gene_id %in% ensembl, ]

# Filter and transpose the data matrix
data <- assay(df)
data <- data[rownames(data) %in% ensembl, colnames(data) %in% rownames(colData)]
rownames(data) <- rowData$external_gene_name
data <- t(data)

# Create DESeq2 dataset
dds <- DESeqDataSetFromMatrix(countData = t(data), colData = colData, design = ~1)

# Normalize the data using DESeq2
dds <- estimateSizeFactors(dds)
normalized_counts <- counts(dds, normalized = TRUE)

# Transpose back the normalized counts
normalized_data <- t(normalized_counts)

# Combine normalized gene expression data with clinical data
data <- cbind(normalized_data, colData)

# Remove columns with all zero values
non_zero_cols <- apply(data, 2, function(col) any(col != 0))
data <- data[, non_zero_cols]
write.csv(data, "expression_data.csv", row.names = FALSE)

# Calculate hazard ratios and p-values for each gene
results <- data.frame(
  gene = colnames(data)[1:(ncol(data) - 2)],
  median = numeric(ncol(data) - 2),
  hazard_ratio = numeric(ncol(data) - 2),
  p_value = numeric(ncol(data) - 2),
  coefficient = numeric(ncol(data) - 2)
)

for (i in seq_len(ncol(data) - 2)) {
  gene_data <- data[, i]
  med <- median(gene_data)
  gene_data_binarized <- ifelse(gene_data > med, 1, 0)
  cox.mod <- coxph(Surv(survival_time, vital_status) ~ gene_data_binarized, data = data)
  summary_cox <- summary(cox.mod)
  
  results$median[i] <- med
  results$hazard_ratio[i] <- exp(cox.mod$coefficients)
  results$p_value[i] <- summary_cox$coefficients[1, "Pr(>|z|)"]
  results$coefficient[i] <- cox.mod$coefficients
}

# Filter results based on p-value < 0.05
filtered_results <- results[results$p_value < 0.05, ]

# Classify genes as good or bad prognostic markers
filtered_results <- filtered_results %>%
  mutate(marker_type = ifelse(hazard_ratio > 1, "BPM", "GPM"))

# Sort by hazard ratio
good_prognostic_markers <- filtered_results %>% filter(marker_type == "GPM") %>% arrange(hazard_ratio) %>% head(10)
bad_prognostic_markers <- filtered_results %>% filter(marker_type == "BPM") %>% arrange(desc(hazard_ratio)) %>% head(10)

write.csv(results, "results.csv", row.names = FALSE)
write.csv(good_prognostic_markers, "good_prognostic_markers.csv", row.names = FALSE)
write.csv(bad_prognostic_markers, "bad_prognostic_markers.csv", row.names = FALSE)

# Create a copy of the data for risk classification
high_low_risk <- data[, 1:(ncol(data) - 2)]

# Vectorized operation to classify risk
for (j in 1:ncol(high_low_risk)) {
  high_low_risk[, j] <- ifelse(data[, j] > median(data[, j]), "High", "Low")
}

High_count <- rowSums(high_low_risk == "High")
Low_count <- rowSums(high_low_risk == "Low")

# Create a new column for risk classification in data
data <- cbind(data, High_count, Low_count)
data$risk <- ifelse(data$High_count > ncol(high_low_risk) / 2, "High_risk", "Low_risk")

# Calculate the C-index using Hmisc package
cox.mod.final <- coxph(Surv(survival_time, vital_status) ~ risk, data = data)
predicted_survival <- predict(cox.mod.final, newdata = data, type = "risk")
c_index <- rcorr.cens(predicted_survival, Surv(data$survival_time, data$vital_status))["C Index"]
print(paste("C-index: ", c_index))

# Plot Kaplan-Meier Curves and perform log-rank test
fit <- survfit(Surv(survival_time, vital_status) ~ risk, data = data)
ggsurvplot(fit, data = data, pval = TRUE, risk.table = TRUE, 
           title = "Kaplan-Meier Survival Curves",
           legend.title = "Risk Group",
           legend.labs = levels(data$risk))

# Log-rank test
log_rank_test <- survdiff(Surv(survival_time, vital_status) ~ risk, data = data)
log_rank_p_value <- 1 - pchisq(log_rank_test$chisq, length(log_rank_test$n) - 1)
print(paste("Log-rank test p-value: ", log_rank_p_value))

# Save expression data for good and bad prognostic markers
save_marker_data <- function(marker_data, marker_names, filename) {
  exp_data <- data[, colnames(data) %in% marker_names]
  exp_data <- cbind(data[, c("vital_status", "survival_time", "risk")], exp_data)
  write.csv(exp_data, filename, row.names = FALSE)
}

save_marker_data(data, good_prognostic_markers$gene, "GPM_exp_data.csv")
save_marker_data(data, bad_prognostic_markers$gene, "BPM_exp_data.csv")
