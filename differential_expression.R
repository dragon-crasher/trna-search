library(DESeq2)

# Load count data
countData <- read.csv("/mnt/d/bioinformatics/RNAseq_pipeline/diffexp/merged_tRF_data.csv", check.names = FALSE)

# Remove the first unnamed column
countData <- countData[, -1]

# Set tRF sequence as row names
rownames(countData) <- countData[, "tRF sequence"]

# Remove the tRF sequence column
countData <- countData[, !names(countData) %in% "tRF sequence"]

# Convert to matrix
countdata <- as.matrix(countData)

# Load metadata
coldata <- read.csv("/mnt/d/bioinformatics/RNAseq_pipeline/diffexp/colData.csv", row.names = 1, check.names = FALSE)

# Find common samples and subset
common_samples <- intersect(colnames(countdata), rownames(coldata))
countdata <- countdata[, common_samples]
coldata <- coldata[common_samples, , drop = FALSE]  # drop = FALSE to preserve colData as a data frame

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
dds <- DESeq(dds)

# Extract results
res <- results(dds)

# Convert results to a data frame
res_df <- as.data.frame(res)

# Sort by log2FoldChange
res_df <- res_df[order(res_df$log2FoldChange, decreasing = TRUE), ]

# Export to CSV
write.csv(res_df, file = "/mnt/d/bioinformatics/RNAseq_pipeline/diffexp/DE_results.csv")
