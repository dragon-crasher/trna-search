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

# Set output directory - make it configurable
output_dir <- if (length(args) >= 3) args[3] else getwd()
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
  cat("Created output directory:", output_dir, "\n")
}

cat("Reading count data from:", count_data_file, "\n")
cat("Reading colData from:", coldata_file, "\n")
cat("Output directory:", output_dir, "\n")

# Load count data
countData <- read.csv(count_data_file, check.names = FALSE)

# Remove the first unnamed column if present
if (names(countData)[1] == "") {
  countData <- countData[, -1]
}

# Set 'License Plate' as row names
if (!"License Plate" %in% colnames(countData)) {
  stop("Column 'License Plate' not found in count data.")
}
gene_names <- countData$'License Plate'  # Store gene names before conversion
rownames(countData) <- gene_names

# Remove identifier columns not part of counts
countData <- countData[, !(names(countData) %in% c("License Plate", "tRF sequence", "tRF type(s)"))]

# Convert to matrix and ensure numeric
countdata <- as.matrix(countData)
countdata <- apply(countdata, 2, function(x) {
  x_num <- as.numeric(x)
  if (any(is.na(x_num) & !is.na(x))) {
    stop("Non-numeric values found in count data that cannot be converted to numbers")
  }
  return(x_num)
})

# Restore rownames after apply
rownames(countdata) <- gene_names

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

# Check for NA in condition column
if (any(is.na(coldata$condition))) {
  stop("NA values found in 'condition' column of colData. Please fix before running.")
}

# Ensure condition column is a factor with consistent ordering
# Sort conditions alphabetically for consistent results
condition_levels <- sort(unique(coldata$condition))
coldata$condition <- factor(coldata$condition, levels = condition_levels)

# Set reference level - try common reference names, otherwise use first level
reference_candidates <- c("control", "Control", "CONTROL", "WT", "wt", "normal", "Normal")
reference_level <- NULL
for (ref in reference_candidates) {
  if (ref %in% levels(coldata$condition)) {
    reference_level <- ref
    break
  }
}

if (!is.null(reference_level)) {
  coldata$condition <- relevel(coldata$condition, ref = reference_level)
  cat("Set reference level to:", reference_level, "\n")
} else {
  cat("No standard reference level found. Using first level:", levels(coldata$condition)[1], "\n")
}

# Trim whitespace from sample names in both datasets
colnames(countdata) <- trimws(colnames(countdata))
rownames(coldata) <- trimws(rownames(coldata))

# Find common samples and subset
common_samples <- intersect(colnames(countdata), rownames(coldata))
if (length(common_samples) == 0) {
  stop("No common samples found between count data and colData.")
}

cat("Found", length(common_samples), "common samples out of", 
    ncol(countdata), "count samples and", nrow(coldata), "metadata samples\n")

countdata <- countdata[, common_samples, drop = FALSE]
coldata <- coldata[common_samples, , drop = FALSE]

# Replace NA with zero
countdata[is.na(countdata)] <- 0

# Remove zero-count samples and genes
zero_samples <- colnames(countdata)[colSums(countdata) == 0]
if (length(zero_samples) > 0) {
  cat("Warning: Samples with zero total counts will be removed:\n")
  print(zero_samples)
  countdata <- countdata[, colSums(countdata) > 0, drop = FALSE]
  coldata <- coldata[colnames(countdata), , drop = FALSE]
}

zero_genes <- rownames(countdata)[rowSums(countdata) == 0]
if (length(zero_genes) > 0) {
  cat("Warning:", length(zero_genes), "genes with zero total counts will be removed\n")
  countdata <- countdata[rowSums(countdata) > 0, , drop = FALSE]
}

# Final checks before DESeq2
if (sum(is.na(countdata)) > 0) {
  stop("NA values found in countdata after cleaning. Please fix.")
}

if (!all(colnames(countdata) == rownames(coldata))) {
  stop("Sample names in countdata and coldata do not match exactly after cleaning.")
}

cat("Final dimensions - countdata:", dim(countdata), "coldata:", dim(coldata), "\n")
cat("Conditions found:", paste(levels(coldata$condition), collapse = ", "), "\n")

# Create DESeq2 object
dds <- DESeqDataSetFromMatrix(countData = countdata,
                              colData = coldata,
                              design = ~ condition)

# Filter out genes with low counts (sum <= 1 across all samples)
keep <- rowSums(counts(dds)) > 1
dds <- dds[keep, ]
cat("After filtering low counts:", nrow(dds), "genes retained\n")

# Run DESeq2 analysis
cat("Running DESeq2 analysis...\n")
dds <- estimateSizeFactors(dds, type = 'poscounts')
dds <- DESeq(dds)

# Check available result names
available_results <- resultsNames(dds)
cat("Available results from DESeq2:\n")
print(available_results)

# Get all unique conditions
conditions <- levels(coldata$condition)
cat("Performing pairwise comparisons for", length(conditions), "conditions\n")

# Generate all pairwise combinations of conditions
pairs <- combn(conditions, 2, simplify = FALSE)

# Initialize vector to collect saved file paths
saved_files <- character()

# Loop through each pair and perform DE analysis
for (i in seq_along(pairs)) {
  pair <- pairs[[i]]
  cat("Running DE analysis", i, "of", length(pairs), ":", pair[2], "vs", pair[1], "\n")
  
  tryCatch({
    # Use contrast method which is more reliable
    res <- results(dds, contrast = c("condition", pair[2], pair[1]), alpha = 0.05)
    
    # For shrinkage, try to find the appropriate coefficient name
    # DESeq2 creates names like "condition_B_vs_A" where A is reference
    possible_coef_names <- c(
      paste0("condition_", pair[2], "_vs_", pair[1]),
      paste0("condition", pair[2], "vs", pair[1]),
      paste0("condition.", pair[2], ".vs.", pair[1])
    )
    
    coef_name <- NULL
    for (name in possible_coef_names) {
      if (name %in% available_results) {
        coef_name <- name
        break
      }
    }
    
    # Apply shrinkage if coefficient found, otherwise use original results
    if (!is.null(coef_name)) {
      cat("Using coefficient:", coef_name, "for shrinkage\n")
      res_shrunk <- lfcShrink(dds, coef = coef_name, type = "apeglm")
    } else {
      cat("Could not find matching coefficient for shrinkage, using original results\n")
      res_shrunk <- res
    }
    
    # Convert to dataframe and add linear fold change
    res_df <- as.data.frame(res_shrunk)
    res_df$linearFoldChange <- 2^res_df$log2FoldChange
    res_df <- res_df[order(abs(res_df$log2FoldChange), decreasing = TRUE), ]
    
    # Create safe filename
    safe_pair1 <- gsub("[^a-zA-Z0-9]", "_", pair[1])
    safe_pair2 <- gsub("[^a-zA-Z0-9]", "_", pair[2])
    output_file <- file.path(output_dir, paste0("DE_results_", safe_pair2, "_vs_", safe_pair1, "_", coldata_basename, ".csv"))
    
    # Write results
    write.csv(res_df, file = output_file)
    cat("Saved results to:", output_file, "\n")
    
    # Append to saved files list
    saved_files <- c(saved_files, output_file)
    
  }, error = function(e) {
    cat("Error in analysis for", pair[2], "vs", pair[1], ":", e$message, "\n")
  })
}

# Save session info for reproducibility
sessioninfo_file <- file.path(output_dir, paste0("sessionInfo_", coldata_basename, ".txt"))
tryCatch({
  writeLines(capture.output(sessionInfo()), sessioninfo_file)
  cat("Session info saved to:", sessioninfo_file, "\n")
}, error = function(e) {
  cat("Could not save session info:", e$message, "\n")
})

# Print summary
cat("\n### Analysis Summary ###\n")
cat("Total comparisons attempted:", length(pairs), "\n")
cat("Successful comparisons:", length(saved_files), "\n")
cat("Output files created:\n")
cat(paste(saved_files, collapse = "\n"), "\n")

cat("\nAll pairwise DE analyses completed.\n")