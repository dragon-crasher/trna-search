coldata <- read.csv("/mnt/d/bioinformatics/RNAseq_pipeline/diffexp/colData.csv", row.names = 1)
print(head(colData))
countdata <- as.matrix(read.csv("/mnt/d/bioinformatics/RNAseq_pipeline/diffexp/merged_tRF_data.csv", row.names = 1))

print(dim(countData))
print(head(countData))

all(rownames(colData) == colnames(countData))  # Should return TRUE


# Convert count data to numeric matrix
countData <- countData[, sapply(countData, is.numeric)]
countData <- as.matrix(countData)
