#!/bin/bash
set -euo pipefail

# Check for args and set defaults
if [ $# -lt 3 ]; then
    echo "Usage: $0 <input_fastq_R1> <input_fastq_R2> <parameter_file> [output_dir]"
    echo "Example: $0 sample_R1.fastq sample_R2.fastq parameters.txt sample_paired_results"
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

# Check input files
for f in "$INPUT_FILE_R1" "$INPUT_FILE_R2" "$PARAMETER_FILE" adapters.txt; do
    if [[ ! -f "$f" ]]; then
        echo "ERROR: $f not found! Exiting..."
        exit 1
    fi
done

# Source parameters file
source "$PARAMETER_FILE"

echo "Processing paired files:"
echo "  R1: $INPUT_FILE_R1"
echo "  R2: $INPUT_FILE_R2"
echo "Output directory: $OUTPUT_DIR"
echo "Current shell: $SHELL"

# Import adapters
while IFS='=' read -r name seq; do
    [[ -z "$name" || "$name" == \#* ]] && continue
    name=$(echo "$name" | xargs)
    seq=$(echo "$seq" | xargs)
    seq=${seq#\"}
    seq=${seq%\"}
    declare "$name=$seq"
done < adapters.txt

# Source parameters again in case adapters.txt defines variables used in params
source "$PARAMETER_FILE"

MINT_OUTPUT_DIR="$MAINWORKDIR/MINT/outputs/$OUTPUT_DIR/"
FASTQC_OUTPUT_DIR="$MAINWORKDIR/RNAseq_pipeline/data/fastqc/$OUTPUT_DIR/"

mkdir -p "$MINT_OUTPUT_DIR"
mkdir -p "$FASTQC_OUTPUT_DIR"

# FastQC on raw files
echo "Running FastQC on raw paired files..."
fastqc "$INPUT_FILE_R1" "$INPUT_FILE_R2" -o "$FASTQC_OUTPUT_DIR" -t 20

# Adapter args
ADAPTER_ARGS=""
for adapter_name in $ADAPTER_LIST; do
    adapter_seq="${!adapter_name}"
    ADAPTER_ARGS+=" -a $adapter_seq"
done

ACCESSION=$(basename "$INPUT_FILE_R1" _R1.fastq)
TRIMMED_R1="${ACCESSION}_trimmed_R1.fastq"
TRIMMED_R2="${ACCESSION}_trimmed_R2.fastq"

# Cutadapt paired-end trimming
echo "Running cutadapt for paired-end adapter trimming..."
cutadapt $ADAPTER_ARGS -m "$MINLEN" -M "$MAXLEN" -q "$MIN_QUALITY" -O 5 -n 1 -j 20 --match-read-wildcards \
    -o "$MAINWORKDIR/SRA/fastq/$ACCESSION/$TRIMMED_R1" \
    -p "$MAINWORKDIR/SRA/fastq/$ACCESSION/$TRIMMED_R2" \
    "$INPUT_FILE_R1" "$INPUT_FILE_R2"

echo "Cutadapt finished running!"

# FastQC on trimmed files
echo "Running FastQC on trimmed paired files..."
fastqc "$MAINWORKDIR/SRA/fastq/$ACCESSION/$TRIMMED_R1" "$MAINWORKDIR/SRA/fastq/$ACCESSION/$TRIMMED_R2" -o "$FASTQC_OUTPUT_DIR" -t 20

# Run MINTmap on trimmed R1
echo "Running MINTmap on trimmed R1 file..."
MINT_R1_OUT="$MINT_OUTPUT_DIR/R1"
mkdir -p "$MINT_R1_OUT"
cd "$MAINWORKDIR/MINT"
./MINTmap.pl -f "$MAINWORKDIR/SRA/fastq/$ACCESSION/$TRIMMED_R1" -p "$MINT_R1_OUT"

# Run MINTmap on trimmed R2
echo "Running MINTmap on trimmed R2 file..."
MINT_R2_OUT="$MINT_OUTPUT_DIR/R2"
mkdir -p "$MINT_R2_OUT"
./MINTmap.pl -f "$MAINWORKDIR/SRA/fastq/$ACCESSION/$TRIMMED_R2" -p "$MINT_R2_OUT"

# Merge tRF counts from R1 and R2
echo "Merging MINTmap tRF counts from R1 and R2..."


# Python script to merge counts by summing
echo "Pipeline completed successfully!"
echo "Results stored in:"
echo "  FastQC: $FASTQC_OUTPUT_DIR"
echo "  MINTmap R1 output: $MINT_R1_OUT"
echo "  MINTmap R2 output: $MINT_R2_OUT"
echo "  Merged counts: $MERGED_COUNTS"
