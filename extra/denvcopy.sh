#!/bin/bash
set -e

# Configuration
INPUT_DIR="/mnt/d/bioinformatics/RNAseq_pipeline/data/raw"
OUTPUT_DIR="/mnt/d/bioinformatics/RNAseq_pipeline/results"
FASTQC_DIR="$OUTPUT_DIR/fastqc"
TRIMMED_DIR="$OUTPUT_DIR/trimmed"
MINTMAP_DIR="$OUTPUT_DIR/mintmap"
DESEQ_DIR="$OUTPUT_DIR/deseq2"
THREADS=4

# Create output directories
mkdir -p "$FASTQC_DIR" "$TRIMMED_DIR" "$MINTMAP_DIR" "$DESEQ_DIR"

# Log file
LOG_FILE="$OUTPUT_DIR/pipeline_log.txt"
echo "Starting tRNA fragment analysis pipeline at $(date)" > "$LOG_FILE"

# Function to process paired-end sample
process_sample() {
    local SAMPLE_NAME=$1
    local R1="$INPUT_DIR/${SAMPLE_NAME}_R1.fastq.gz"
    local R2="$INPUT_DIR/${SAMPLE_NAME}_R2.fastq.gz"
    
    echo "Processing sample: $SAMPLE_NAME" | tee -a "$LOG_FILE"
    
    # Step 1: Quality assessment with FastQC
    echo "Running FastQC on $SAMPLE_NAME..." | tee -a "$LOG_FILE"
    fastqc -o "$FASTQC_DIR" -t "$THREADS" "$R1" "$R2"
    
    # Step 2: Adapter trimming with Cutadapt v4.0
    echo "Running Cutadapt on $SAMPLE_NAME..." | tee -a "$LOG_FILE"
    local TRIMMED_R1="$TRIMMED_DIR/${SAMPLE_NAME}_R1_trimmed.fastq.gz"
    local TRIMMED_R2="$TRIMMED_DIR/${SAMPLE_NAME}_R2_trimmed.fastq.gz"
    
    # Adapter sequences - replace with the actual adapters used in your study
    local ADAPTER_R1="AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"
    local ADAPTER_R2="AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"
    
    cutadapt \
        -a "$ADAPTER_R1" \
        -A "$ADAPTER_R2" \
        --minimum-length 15 \
        --maximum-length 50 \
        -o "$TRIMMED_R1" \
        -p "$TRIMMED_R2" \
        --cores="$THREADS" \
        "$R1" "$R2" \
        > "$TRIMMED_DIR/${SAMPLE_NAME}_cutadapt_report.txt"
    
    # Step 3: Run MINTmap for tRNA fragment identification and quantification
    echo "Running MINTmap on $SAMPLE_NAME..." | tee -a "$LOG_FILE"
    local MINTMAP_OUTPUT="$MINTMAP_DIR/${SAMPLE_NAME}"
    mkdir -p "$MINTMAP_OUTPUT"
    
    # Assuming MINTmap is installed and in PATH
    # Note: You may need to adjust this command based on MINTmap's actual interface
    mintmap \
        --input "$TRIMMED_R1,$TRIMMED_R2" \
        --output "$MINTMAP_OUTPUT" \
        --paired \
        --threads "$THREADS"
    
    echo "Completed processing for $SAMPLE_NAME" | tee -a "$LOG_FILE"
}

# Process all samples
for SAMPLE_PATH in "$INPUT_DIR"/*_R1.fastq.gz; do
    SAMPLE_NAME=$(basename "$SAMPLE_PATH" | sed 's/_R1.fastq.gz//')
    process_sample "$SAMPLE_NAME"
done

# Step 4: Differential Expression Analysis with DESeq2
# This would typically be done in R - creating an R script
echo "Running differential expression analysis with DESeq2..." | tee -a "$LOG_FILE"

cat > "$DESEQ_DIR/run_deseq2.R" << 'EOF'
#!/usr/bin/env Rscript

# Load required libraries
library(DESeq2)
library(tidyverse)

# Set working directory to the MINTmap output
mintmap_dir <- commandArgs(trailingOnly = TRUE)[1]
output_dir <- commandArgs(trailingOnly = TRUE)[2]

# Read sample information
# Create a sample information file manually or modify this part accordingly
samples <- data.frame(
  sample = list.files(mintmap_dir, pattern = ".*"),
  condition = factor(c("control", "control", "treated", "treated"))  # Modify based on your experimental design
)

# Create count matrix from MINTmap output
count_files <- file.path(mintmap_dir, samples$sample, "fragment_counts.txt")
names(count_files) <- samples$sample

# Function to read MINTmap count files
read_mintmap <- function(file) {
  counts <- read.table(file, header = TRUE, sep = "\t")
  return(counts$Count)  # Adjust column name if needed
}

# Create count matrix
counts_data <- do.call(cbind, lapply(count_files, read_mintmap))
rownames(counts_data) <- read.table(count_files[1], header = TRUE, sep = "\t")$Fragment  # Adjust column name if needed

# Create DESeq2 dataset
dds <- DESeqDataSetFromMatrix(
  countData = counts_data,
  colData = samples,
  design = ~ condition
)

# Run DESeq2 analysis
dds <- DESeq(dds)

# Get results
res <- results(dds, contrast = c("condition", "treated", "control"))

# Filter results by adjusted p-value and fold change
significant_results <- res[!is.na(res$padj) & res$padj < 0.05 & abs(res$log2FoldChange) >= 2, ]

# Save results
write.csv(as.data.frame(res), file = file.path(output_dir, "all_results.csv"))
write.csv(as.data.frame(significant_results), file = file.path(output_dir, "significant_results.csv"))

# Create a basic volcano plot
pdf(file.path(output_dir, "volcano_plot.pdf"), width = 10, height = 8)
plot(res$log2FoldChange, -log10(res$padj), 
     pch = 20, col = ifelse(res$padj < 0.05 & abs(res$log2FoldChange) >= 2, "red", "black"),
     xlab = "Log2 Fold Change", ylab = "-Log10 Adjusted P-value",
     main = "Volcano Plot of Differentially Expressed tRNA Fragments")
abline(v = c(-2, 2), lty = 2)
abline(h = -log10(0.05), lty = 2)
dev.off()

# Print summary
cat("DESeq2 analysis completed.\n")
cat("Total number of tRNA fragments analyzed:", nrow(res), "\n")
cat("Number of significantly differentially expressed fragments:", nrow(significant_results), "\n")
EOF

# Run the R script
Rscript "$DESEQ_DIR/run_deseq2.R" "$MINTMAP_DIR" "$DESEQ_DIR"

echo "Pipeline completed at $(date)" | tee -a "$LOG_FILE"
echo "Results are available in $OUTPUT_DIR"