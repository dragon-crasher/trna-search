library(DESeq2)

# Load count data and metadata
countdata <- as.matrix(read.csv("/mnt/d/bioinformatics/RNAseq_pipeline/diffexp/merged_tRF_data.csv", row.names = 1))
storage.mode(countdata) <- "numeric"

coldata <- read.csv("/mnt/d/bioinformatics/RNAseq_pipeline/diffexp/colData.csv", row.names = 1)

# Ensure colData and countData have matching sample names in the correct order
common_samples <- intersect(colnames(countdata), rownames(coldata))
countdata <- countdata[, common_samples]
coldata <- coldata[common_samples, ]

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

# Extract results and export to CSV
res <- results(dds)
write.csv(as.data.frame(res), file = "/mnt/d/bioinformatics/RNAseq_pipeline/diffexp/DE_results.csv")
