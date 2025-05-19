library(DESeq2)

# Load count data
countData <- read.csv("/mnt/d/bioinformatics/RNAseq_pipeline/data/CLLmerged_mint_files.csv", check.names = FALSE)

# Remove the first unnamed column
if (names(countData)[1] == "") {
  countData <- countData[, -1]
}
# Set tRF sequence as row names
rownames(countData) <- countData$'MINTbase Unique ID'

# Remove the tRF sequence column
countData <- countData[, !(names(countData) %in% c("MINTbase Unique ID", "tRF sequence", "tRF type(s)"))]

# Convert to matrix
countdata <- as.matrix(countData)

# Load metadata
coldata <- read.csv("/mnt/d/bioinformatics/RNAseq_pipeline/diffexp/colData.csv", row.names = 1, check.names = FALSE)

# Clean up condition column: replace spaces and special characters with underscores
coldata$condition <- gsub(" ", "_", coldata$condition)
coldata$condition <- gsub("[^a-zA-Z0-9_.]", "_", coldata$condition)

# Ensure condition column is a factor
coldata$condition <- factor(coldata$condition)

# Find common samples and subset
common_samples <- intersect(colnames(countdata), rownames(coldata))
countdata <- countdata[, common_samples]
coldata <- coldata[common_samples, , drop = FALSE]  # drop = FALSE to preserve colData as a data frame
countdata[is.na(countdata)] <- 0

cat("Samples in countdata:\n")
print(colnames(countdata))

cat("\nSamples in coldata:\n")
print(rownames(coldata))


# Ensure countdata and coldata are not empty after subsetting
if (ncol(countdata) == 0 || nrow(coldata) == 0) {
    stop("Error: countdata or coldata is empty after subsetting. Check sample names.")
}

# Verify dimensions
print(dim(countdata))
print(dim(coldata))

# Verify that colData and countData have matching sample names
print(identical(colnames(countdata), rownames(coldata)))

# Create DESeq2 object
dds <- DESeqDataSetFromMatrix(countData = countdata,
                                colData = coldata,
                                design = ~ condition)

# Run DESeq2 analysis
dds <- estimateSizeFactors(dds, type='poscounts')
dds <- DESeq(dds)
# Extract results for healthy vs indolent comparison
res_indolent <- results(dds, contrast = c("condition", "INDOLENT_CLL", "HEALTHY"))

# Extract results for healthy vs aggressive comparison
res_aggressive <- results(dds, contrast = c("condition", "AGGRESSIVE_CLL", "HEALTHY"))

# Convert results to data frames
res_indolent_df <- as.data.frame(res_indolent)
res_aggressive_df <- as.data.frame(res_aggressive)

# Calculate linear fold change (2^log2FoldChange)
res_indolent_df$linearFoldChange <- 2^res_indolent_df$log2FoldChange
res_aggressive_df$linearFoldChange <- 2^res_aggressive_df$log2FoldChange

# Sort by linear fold change (absolute value)
res_indolent_df <- res_indolent_df[order(abs(res_indolent_df$linearFoldChange), decreasing = TRUE), ]
res_aggressive_df <- res_aggressive_df[order(abs(res_aggressive_df$linearFoldChange), decreasing = TRUE), ]

# Export to CSV
write.csv(res_indolent_df, file = "/mnt/d/bioinformatics/RNAseq_pipeline/data/DE_results_healthy_vs_indolent.csv")
write.csv(res_aggressive_df, file = "/mnt/d/bioinformatics/RNAseq_pipeline/data/DE_results_healthy_vs_aggressive.csv")
