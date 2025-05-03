#!/bin/bash

# Check args and set defaults
if [ $# -lt 3 ]; then
    echo "Usage: $0 <input_fastq_R1> <input_fastq_R2> <parameter_file> [output_dir]"
    exit 1
fi

INPUT_FILE_R1=$1
INPUT_FILE_R2=$2
PARAMETER_FILE=$3

if [ $# -ge 4 ]; then
    OUTPUT_DIR=$4
else
    OUTPUT_DIR=$(basename "$INPUT_FILE_R1" _R1.fastq)
fi

if [[ ! -f "$INPUT_FILE_R1" ]]; then
    echo "ERROR: $INPUT_FILE_R1 not found! Exiting..."
    exit 1
fi

if [[ ! -f "$INPUT_FILE_R2" ]]; then
    echo "ERROR: $INPUT_FILE_R2 not found! Exiting..."
    exit 1
fi

if [[ ! -f "$PARAMETER_FILE" ]]; then
    echo "ERROR: $PARAMETER_FILE not found! Exiting..."
    exit 1
fi

if [[ ! -f "adapters.txt" ]]; then
    echo "ERROR: adapters.txt not found! Exiting..."
    exit 1
fi

source adapters.txt
source "$PARAMETER_FILE"

FASTQC_OUTPUT_DIR="/mnt/d/bioinformatics/RNAseq_pipeline/data/fastqc/$OUTPUT_DIR"
mkdir -p "$FASTQC_OUTPUT_DIR"

echo "Running FastQC on raw files..."
fastqc "$INPUT_FILE_R1" "$INPUT_FILE_R2" -o "$FASTQC_OUTPUT_DIR" -t 6

# Prepare adapter args as before
ADAPTER_ARGS=""
for adapter_name in $ADAPTER_LIST; do
    adapter_seq="${!adapter_name}"
    ADAPTER_ARGS+=" -a $adapter_seq"
done

echo "Running cutadapt for paired-end trimming..."
cutadapt $ADAPTER_ARGS -m "$MINLEN" -M "$MAXLEN" -q "$MIN_QUALITY" -O 5 --match-read-wildcards \
    -o trimmed_R1.fastq -p trimmed_R2.fastq "$INPUT_FILE_R1" "$INPUT_FILE_R2"

echo "Cutadapt finished!"

echo "Running FastQC on trimmed files..."
fastqc trimmed_R1.fastq trimmed_R2.fastq -o "$FASTQC_OUTPUT_DIR" -t 6

# For MINTmap, check if paired-end supported; if not, process trimmed_R1.fastq or both separately
# Example (if single-end only):
MINT_OUTPUT_DIR="/mnt/d/bioinformatics/MINT/outputs/$OUTPUT_DIR"
mkdir -p "$MINT_OUTPUT_DIR"

echo "Running MINTmap on trimmed_R1.fastq..."
cd /mnt/d/bioinformatics/MINT
./MINTmap.pl -f trimmed_R1.fastq -p "$MINT_OUTPUT_DIR"

# Report runtime etc.
duration=$SECONDS
echo "Pipeline completed in $(($duration/60)) minutes and $(($duration % 60)) seconds."
echo "Results stored in:"
echo "  FastQC: $FASTQC_OUTPUT_DIR"
echo "  MINTmap: $MINT_OUTPUT_DIR"