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
library(survMisc)

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


results <- data.frame(Gene = character(), HR = numeric(), p_value = numeric(), CI_lower = numeric(), CI_upper = numeric())

for (gene in colnames(data)[1:(ncol(data) - 2)]) {
  med <- median(data[[gene]])
  gene_data_binarized <- ifelse(data[[gene]] > med, 1, 0)
  cox_model <- coxph(Surv(survival_time, vital_status) ~ gene_data_binarized, data = data)
  cox_summary <- summary(cox_model)
  results <- rbind(results, data.frame(
    Gene = gene,
    HR = exp(cox_model$coefficients),
    p_value = cox_summary$coefficients[,"Pr(>|z|)"],
    CI_lower = exp(confint(cox_model))[1],
    CI_upper = exp(confint(cox_model))[2]
  ))
}

# Filter significant genes
significant_genes <- results[results$p_value < 0.05, ]

# Calculate Prognostic Index (PI)
coefficients <- significant_genes$HR
names(coefficients) <- significant_genes$Gene
data$PI <- rowSums(sapply(names(coefficients), function(gene) data[[gene]] * coefficients[gene]))

# Ensure survival_time and vital_status are numeric
data$survival_time <- as.numeric(data$survival_time)
data$vital_status <- as.numeric(data$vital_status)

# Determine the optimal cutoff for PI
cutr <- cutp(coxph(Surv(survival_time, vital_status) ~ PI, data = data))[[1]][]
cutr<-cutr[,1]
# Assign risk groups based on the cutoff
data$risk <- ifelse(data$PI >= cutr, "High_risk", "Low_risk")

# Fit a survival curve using the Kaplan-Meier estimator for the risk groups
km_fit <- survfit(Surv(survival_time, vital_status) ~ risk, data = data)

# Plot the Kaplan-Meier survival curves
ggsurvplot(km_fit, data = data, pval = TRUE, risk.table = TRUE)

# Evaluate the model performance using Hazard Ratios and Concordance index
cox_model_pi <- coxph(Surv(survival_time, vital_status) ~ PI, data = data)
summary_cox_pi <- summary(cox_model_pi)

# Calculate Concordance index (C-index)
C_index <- summary_cox_pi$concordance[1]

# Print model evaluation metrics
print(paste("Hazard Ratio:", summary_cox_pi$coefficients[1, "exp(coef)"]))
print(paste("95% CI:", summary_cox_pi$conf.int[1, "lower .95"], "-", summary_cox_pi$conf.int[1, "upper .95"]))
print(paste("p-value:", summary_cox_pi$coefficients[1, "Pr(>|z|)"]))
print(paste("Concordance Index (C):", C_index))

# Log-rank test
log_rank_test <- survdiff(Surv(survival_time, vital_status) ~ risk, data = data)
log_rank_p_value <- 1 - pchisq(log_rank_test$chisq, length(log_rank_test$n) - 1)
print(paste("Log-rank test p-value: ", log_rank_p_value))
write.csv(significant_genes, "significant_genes.csv", row.names = FALSE)
