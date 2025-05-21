#!/usr/bin/env Rscript

library(DESeq2)

# Get command-line arguments (filenames only)
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  stop("Usage: Rscript deseq_analysis_all_pairs.R <count_data_filename> <colData_filename>\n",
       "Example: Rscript deseq_analysis_all_pairs.R CLLmerged_mint_files.csv colData.csv")
}

# Join filenames with directories
count_data_file <- args[1]
coldata_file <- args[2]

# Extract basename of coldata file without extension for output naming
coldata_basename <- tools::file_path_sans_ext(basename(coldata_file))

cat("Reading count data from:", count_data_file, "\n")
cat("Reading colData from:", coldata_file, "\n")

# Load count data
countData <- read.csv(count_data_file, check.names = FALSE)

# Remove the first unnamed column if present
if (names(countData)[1] == "") {
  countData <- countData[, -1]
}

# Set tRF sequence as row names
if (!"License Plate" %in% colnames(countData)) {
  stop("Column 'License Plate' not found in count data.")
}
rownames(countData) <- countData$'License Plate'

# Remove identifier columns not part of counts
countData <- countData[, !(names(countData) %in% c("License Plate", "tRF sequence", "tRF type(s)"))]

# Convert to matrix
countdata <- as.matrix(countData)

# Load metadata
coldata <- read.csv(coldata_file, check.names = FALSE)

# Set rownames from Sample column and remove Sample column
if (!"Sample" %in% colnames(coldata)) {
  stop("Column 'Sample' not found in colData.")
}
rownames(coldata) <- coldata$Sample
coldata$Sample <- NULL

# Remove duplicate samples in colData
if (any(duplicated(rownames(coldata)))) {
  dup_samples <- rownames(coldata)[duplicated(rownames(coldata))]
  cat("Warning: Duplicate samples found in colData and will be removed:\n")
  print(dup_samples)
  coldata <- coldata[!duplicated(rownames(coldata)), ]
}

# Clean up condition column: replace spaces and special characters with underscores
if (!"condition" %in% colnames(coldata)) {
  stop("Column 'condition' not found in colData")
}
coldata$condition <- gsub(" ", "_", coldata$condition)
coldata$condition <- gsub("[^a-zA-Z0-9_.]", "_", coldata$condition)

# Ensure condition column is a factor
coldata$condition <- factor(coldata$condition)

# Trim whitespace from sample names in both datasets
colnames(countdata) <- trimws(colnames(countdata))
rownames(coldata) <- trimws(rownames(coldata))

# Find common samples and subset
common_samples <- intersect(colnames(countdata), rownames(coldata))
if (length(common_samples) == 0) {
  stop("No common samples found between count data and colData.")
}

countdata <- countdata[, common_samples, drop = FALSE]
coldata <- coldata[common_samples, , drop = FALSE]

# Replace NA counts with zero
countdata[is.na(countdata)] <- 0

cat("Samples in countdata:\n")
print(colnames(countdata))

cat("\nSamples in coldata:\n")
print(rownames(coldata))

# Ensure countdata and coldata are not empty after subsetting
if (ncol(countdata) == 0 || nrow(coldata) == 0) {
  stop("Error: countdata or coldata is empty after subsetting. Check sample names.")
}

cat("Dimensions of countdata:", dim(countdata), "\n")
cat("Dimensions of coldata:", dim(coldata), "\n")
cat("Sample names match:", identical(colnames(countdata), rownames(coldata)), "\n")

# Create DESeq2 object
dds <- DESeqDataSetFromMatrix(countData = countdata,
                              colData = coldata,
                              design = ~ condition)

# Filter out genes with low counts (sum <= 1 across all samples)
keep <- rowSums(counts(dds)) > 1
dds <- dds[keep, ]

# Run DESeq2 analysis
dds <- estimateSizeFactors(dds, type = 'poscounts')
dds <- DESeq(dds)

# Get all unique conditions
conditions <- levels(coldata$condition)

# Generate all pairwise combinations of conditions
pairs <- combn(conditions, 2, simplify = FALSE)

# Loop through each pair and perform DE analysis
for (pair in pairs) {
  cat("Running DE analysis for:", pair[1], "vs", pair[2], "\n")
  
  # Run results with contrast (optional, for p-values etc.)
  res <- results(dds, contrast = c("condition", pair[1], pair[2]), alpha = 0.05)
  
  # Construct coef name: assuming pair[1] is reference, pair[2] is compared
  coef_name <- paste0("condition_", pair[2], "_vs_", pair[1])
  
  # Shrink fold changes using coef
  res_shrunk <- lfcShrink(dds, coef = coef_name, type = "apeglm")
  
  res_df <- as.data.frame(res_shrunk)
  res_df$linearFoldChange <- 2^res_df$log2FoldChange
  res_df <- res_df[order(abs(res_df$linearFoldChange), decreasing = TRUE), ]
  
  safe_pair1 <- gsub("[^a-zA-Z0-9]", "_", pair[1])
  safe_pair2 <- gsub("[^a-zA-Z0-9]", "_", pair[2])
  output_file <- paste0("/raid/anirudh/bioinformatics/RNAseq_pipeline/data/DE_results_", safe_pair1, "_vs_", safe_pair2, "_", coldata_basename, ".csv")
  
  write.csv(res_df, file = output_file)
  cat("Saved results to:", output_file, "\n")
}


# Save session info for reproducibility
sessioninfo_file <- paste0("/raid/anirudh/bioinformatics/RNAseq_pipeline/data/sessionInfo_", coldata_basename, ".txt")
writeLines(capture.output(sessionInfo()), sessioninfo_file)
cat("Session info saved to:", sessioninfo_file, "\n")

cat("All pairwise DE analyses completed.\n")
