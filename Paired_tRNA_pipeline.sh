#!/bin/bash
set -euo pipefail

# Check for args and set defaults
if [ $# -lt 4 ]; then
    echo "Usage: $0 <input_fastq_R1> <input_fastq_R2> <parameter_file> <output_dir> <run_cutadapt>"
    echo "Example: $0 sample_R1.fastq sample_R2.fastq parameters.txt sample_paired_results 1"
    exit 1
fi

INPUT_FILE_R1=$1
INPUT_FILE_R2=$2
PARAMETER_FILE=$3
RUN_CUTADAPT=$5

if [ $# -ge 5 ]; then
    OUTPUT_DIR=$4
else
    OUTPUT_DIR=$(basename "$INPUT_FILE_R1" _1.fastq)
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

MINT_OUTPUT_DIR="$MAINWORKDIR/MINT/outputs/$PROJECT_NAME/"
FASTQC_OUTPUT_DIR="$MAINWORKDIR/RNAseq_pipeline/data/fastqc/$OUTPUT_DIR/"

mkdir -p "$FASTQC_OUTPUT_DIR"
mkdir -p "$MINT_OUTPUT_DIR"

# FastQC on raw files
echo "Running FastQC on raw paired files..."
fastqc "$INPUT_FILE_R1" "$INPUT_FILE_R2" -o "$FASTQC_OUTPUT_DIR" -t $cputhreads

# Adapter args
for adapter_name in $ADAPTER3_LIST; do
    if [[ -z "${!adapter_name:-}" ]]; then
        echo "ERROR: Adapter sequence for '$adapter_name' is not set. Check adapters.txt."
        exit 1
    fi
    adapter_seq="${!adapter_name}"
    ADAPTER_ARGS1+=" -a $adapter_seq"

done

for adapter_name in $ADAPTER5_LIST; do
    if [[ -z "${!adapter_name:-}" ]]; then
        echo "ERROR: Adapter sequence for '$adapter_name' is not set. Check adapters.txt."
        exit 1
    fi
    adapter_seq="${!adapter_name}"
    ADAPTER_ARGS1+=" -g $adapter_seq"

done

for adapter_name in $ADAPTER3_LIST2; do
    if [[ -z "${!adapter_name:-}" ]]; then
        echo "ERROR: Adapter sequence for '$adapter_name' is not set. Check adapters.txt."
        exit 1
    fi
    adapter_seq="${!adapter_name}"
    ADAPTER_ARGS2+=" -A $adapter_seq"
done

for adapter_name in $ADAPTER5_LIST2; do
    if [[ -z "${!adapter_name:-}" ]]; then
        echo "ERROR: Adapter sequence for '$adapter_name' is not set. Check adapters.txt."
        exit 1
    fi
    adapter_seq="${!adapter_name}"
    ADAPTER_ARGS2+=" -G $adapter_seq"
done

ACCESSION=$(basename "$INPUT_FILE_R1" _1.fastq)
TRIMMED_R1="${ACCESSION}_trimmed_R1.fastq"
TRIMMED_R2="${ACCESSION}_trimmed_R2.fastq"

if (( RUN_CUTADAPT == 1)); then
    echo "Running cutadapt for paired-end adapter trimming..."
    mkdir -p "$MAINWORKDIR/SRA/fastq/$ACCESSION"

    cutadapt $ADAPTER_ARGS1 $ADAPTER_ARGS2  -m "$MINLEN" -M "$MAXLEN" -q "$MIN_QUALITY" -O 5 -n 1 -j $cpucores --match-read-wildcards \
        -o "$MAINWORKDIR/SRA/fastq/$ACCESSION/$TRIMMED_R1" \
        -p "$MAINWORKDIR/SRA/fastq/$ACCESSION/$TRIMMED_R2" \
        "$INPUT_FILE_R1" "$INPUT_FILE_R2"

    echo "Cutadapt finished running!"

    echo "Running FastQC on trimmed paired files..."
    fastqc "$MAINWORKDIR/SRA/fastq/$ACCESSION/$TRIMMED_R1" "$MAINWORKDIR/SRA/fastq/$ACCESSION/$TRIMMED_R2" -o "$FASTQC_OUTPUT_DIR" -t $cputhreads

    RUNNING_R1="$MAINWORKDIR/SRA/fastq/$ACCESSION/$TRIMMED_R1"
    RUNNING_R2="$MAINWORKDIR/SRA/fastq/$ACCESSION/$TRIMMED_R2"
else
    echo "Skipping cutadapt and using raw input files for downstream analysis."
    RUNNING_R1="$INPUT_FILE_R1"
    RUNNING_R2="$INPUT_FILE_R2"
fi

# Run MINTmap on the files assigned to RUNNING_R1 and RUNNING_R2
echo "Running MINTmap on R1 file..."
MINT_R1_OUT="R1_$MINT_OUTPUT_DIR"

cd "$MINT_OUTPUT_DIR"
#./MINTmap.pl -f "$RUNNING_R1" -p "$MINT_R1_OUT"
conda run -n mintmap38 MINTmap -p "R1_$OUTPUT_DIR" "$RUNNING_R1" 
echo "Running MINTmap on R2 file..."
MINT_R2_OUT="R2_$MINT_OUTPUT_DIR"
#./MINTmap.pl -f "$RUNNING_R2" -p "$MINT_R2_OUT"
conda run -n mintmap38 MINTmap -p "R2_$OUTPUT_DIR" "$RUNNING_R2"


echo "Pipeline completed successfully!"
echo "Results stored in:"
echo "  FastQC: $FASTQC_OUTPUT_DIR"
echo "  MINTmap R1 output: $MINT_R1_OUT"
echo "  MINTmap R2 output: $MINT_R2_OUT"

