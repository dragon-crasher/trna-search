#!/bin/bash

# Accept input dataset filename as an argument
if [ $# -lt 1 ]; then
    echo "Usage: $0 <input_fastq_file> [output_dir]"
    echo "Example: $0 sample1.fastq sample1_results"
    exit 1
fi

INPUT_FILE=$1
FILENAME=$(basename "$INPUT_FILE" .fastq)

# Set output directory, use second argument if provided, otherwise use filename
if [ $# -ge 2 ]; then
    OUTPUT_DIR=$2
else
    OUTPUT_DIR=$FILENAME
fi

echo "Processing file: $INPUT_FILE"
echo "Output directory: $OUTPUT_DIR"
echo "Current shell: $SHELL"

SECONDS=0

# Create output directories
MINT_OUTPUT_DIR="/mnt/d/bioinformatics/MINT/MINTmap-release-v1.0/MINTmap-release-v1.0/outputs/$OUTPUT_DIR"
FASTQC_OUTPUT_DIR="/mnt/d/bioinformatics/RNAseq_pipeline/data/fastqc/$OUTPUT_DIR"
# HISAT2_OUTPUT_DIR="/mnt/d/bioinformatics/RNAseq_pipeline/HISAT2/$OUTPUT_DIR"
# QUANTS_OUTPUT_DIR="/mnt/d/bioinformatics/RNAseq_pipeline/quants/$OUTPUT_DIR"

mkdir -p "$FASTQC_OUTPUT_DIR"
# mkdir -p "$HISAT2_OUTPUT_DIR"
# mkdir -p "$QUANTS_OUTPUT_DIR"
mkdir -p "$MINT_OUTPUT_DIR"

# Change working directory
cd /mnt/d/bioinformatics/RNAseq_pipeline/

# STEP 1: Run fastqc
echo "Running FastQC on input file..."
fastqc "$INPUT_FILE" -o "$FASTQC_OUTPUT_DIR"

# Run trimmomatic to trim reads with poor quality
# echo "Running Trimmomatic..."
# TRIMMED_FILE="/mnt/d/bioinformatics/RNAseq_pipeline/data/${FILENAME}_trimmed.fastq"
# if [[ ! -f "$INPUT_FILE" ]]; then
#     echo "ERROR: Input FASTQ file $INPUT_FILE not found! Exiting..."
#     exit 1
# fi
# java -jar /mnt/d/bioinformatics/Trimmomatic-0.39/trimmomatic-0.39.jar SE -threads 4 "$INPUT_FILE" "$TRIMMED_FILE" TRAILING:10 -phred33

# echo "Trimmomatic finished running!"
#conda activate cutadapt

echo "Running cutadapt..."
TRIMMED_FILE="/mnt/d/bioinformatics/RNAseq_pipeline/data/${FILENAME}_trimmed.fastq"
if [[ ! -f "$INPUT_FILE" ]]; then
    echo "ERROR: Input FASTQ file $INPUT_FILE not found! Exiting..."
    exit 1
fi
cutadapt -q 10 -o "$TRIMMED_FILE" "$INPUT_FILE" --cores=4

echo "cutadapt finished running!"  

# Run FastQC on trimmed file
echo "Running FastQC on trimmed file..."
fastqc "$TRIMMED_FILE" -o "$FASTQC_OUTPUT_DIR"

# # STEP 2: Run HISAT2 (uncommented for full pipeline)
# echo "Running HISAT2 alignment..."
# hisat2 -q --rna-strandness R -x HISAT2/grch38/genome -U "$TRIMMED_FILE" | samtools sort -o "${HISAT2_OUTPUT_DIR}/${FILENAME}_trimmed.bam"

# # STEP 3: Run featureCounts - Quantification
# echo "Running featureCounts..."
# featureCounts -s 2 -a "/mnt/d/bioinformatics/RNAseq_pipeline/hg38/Homo_sapiens.GRCh38.113.chr_patch_hapl_scaff.gtf.gz" -o "${QUANTS_OUTPUT_DIR}/${FILENAME}_featurecounts.txt" "${HISAT2_OUTPUT_DIR}/${FILENAME}_trimmed.bam"

# Run MINTmap
echo "Running MINTmap..."
cd /mnt/d/bioinformatics/MINT/MINTmap-release-v1.0/MINTmap-release-v1.0
./MINTmap.pl -f "$TRIMMED_FILE" -p "$MINT_OUTPUT_DIR"

# Report runtime
duration=$SECONDS
echo "Pipeline completed in $(($duration/60)) minutes and $(($duration % 60)) seconds."
echo "Results stored in:"
echo "  FastQC: $FASTQC_OUTPUT_DIR"
# echo "  HISAT2: $HISAT2_OUTPUT_DIR"
# echo "  featureCounts: $QUANTS_OUTPUT_DIR" 
echo "  MINTmap: $MINT_OUTPUT_DIR"