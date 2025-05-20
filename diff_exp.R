#!/usr/bin/env Rscript

library(DESeq2)

# Fixed directories for input files
count_data_dir <- "/mnt/d/bioinformatics/RNAseq_pipeline/data/"
coldata_dir <- "/mnt/d/bioinformatics/RNAseq_pipeline/diffexp/"

# Get command-line arguments (filenames only)
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  stop("Usage: Rscript deseq_analysis_all_pairs.R <count_data_filename> <colData_filename>\n",
       "Example: Rscript deseq_analysis_all_pairs.R CLLmerged_mint_files.csv colData.csv")
}

# Join filenames with directories
count_data_file <- file.path(count_data_dir, args[1])
coldata_file <- file.path(coldata_dir, args[2])

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
rownames(countData) <- countData$'License Plate'

# Remove identifier columns
countData <- countData[, !(names(countData) %in% c("License Plate", "tRF sequence", "tRF type(s)"))]

# Convert to matrix
countdata <- as.matrix(countData)

# Load metadata
coldata <- read.csv(coldata_file, row.names = 1, check.names = FALSE)

# Clean up condition column: replace spaces and special characters with underscores
coldata$condition <- gsub(" ", "_", coldata$condition)
coldata$condition <- gsub("[^a-zA-Z0-9_.]", "_", coldata$condition)

# Ensure condition column is a factor
coldata$condition <- factor(coldata$condition)

# Find common samples and subset
common_samples <- intersect(colnames(countdata), rownames(coldata))
countdata <- countdata[, common_samples]
coldata <- coldata[common_samples, , drop = FALSE]
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
cat("Dimensions of countdata:", dim(countdata), "\n")
cat("Dimensions of coldata:", dim(coldata), "\n")

# Verify matching sample names
cat("Sample names match:", identical(colnames(countdata), rownames(coldata)), "\n")

# Create DESeq2 object
dds <- DESeqDataSetFromMatrix(countData = countdata,
                              colData = coldata,
                              design = ~ condition)

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
  
  res <- results(dds, contrast = c("condition", pair[1], pair[2]))
  
  # Convert results to data frame
  res_df <- as.data.frame(res)
  
  # Calculate linear fold change
  res_df$linearFoldChange <- 2^res_df$log2FoldChange
  
  # Sort by absolute linear fold change
  res_df <- res_df[order(abs(res_df$linearFoldChange), decreasing = TRUE), ]
  
  # Create safe output filename
  safe_pair1 <- gsub("[^a-zA-Z0-9]", "_", pair[1])
  safe_pair2 <- gsub("[^a-zA-Z0-9]", "_", pair[2])
  output_file <- paste0("/mnt/d/bioinformatics/RNAseq_pipeline/data/DE_results_", safe_pair1, "_vs_", safe_pair2, "_", coldata_basename, ".csv")
  
  # Write results to CSV
  write.csv(res_df, file = output_file)
  
  cat("Saved results to:", output_file, "\n")
}

cat("All pairwise DE analyses completed.\n")
