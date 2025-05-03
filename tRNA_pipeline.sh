#!/bin/bash

#check for args and set defaults
if [ $# -lt 2 ]; then
    echo "Usage: $0 <input_fastq_file> <parameter_file> [output_dir]"
    echo "Example: $0 sample1.fastq parameters.txt sample1_results"
    exit 1
fi
INPUT_FILE=$1
PARAMETER_FILE=$2

if [[ ! -f "$INPUT_FILE" ]]; then
    echo "ERROR: Input FASTQ file $INPUT_FILE not found! Exiting..."
    exit 1
fi

if [[ ! -f "$PARAMETER_FILE" ]]; then
    echo "ERROR: Parameter file $PARAMETER_FILE not found! Exiting..."
    exit 1
fi

if [[ ! -f "adapters.txt" ]]; then
    echo "ERROR: adapters.txt file not found! Exiting..."
    exit 1
fi

# Set output directory, use third argument if provided, otherwise use filename
if [ $# -ge 3 ]; then
    OUTPUT_DIR=$3
else
    OUTPUT_DIR=$(basename "$1" .fastq)
fi

#import list of adapters and parameters
source adapters.txt
source "$PARAMETER_FILE"

echo "Processing file: $INPUT_FILE"
echo "Output directory: $OUTPUT_DIR"
echo "Current shell: $SHELL"

SECONDS=0
# Create output directories
MINT_OUTPUT_DIR="$MAINWORKDIR/MINT/outputs/$OUTPUT_DIR"
FASTQC_OUTPUT_DIR="$MAINWORKDIR/RNAseq_pipeline/data/fastqc/$OUTPUT_DIR"

#mkdir -p "$MINT_OUTPUT_DIR"

#perform quality check
mkdir -p "$FASTQC_OUTPUT_DIR"
echo "Running FastQC on input file..."
first_header=$(zcat "$FILE" | head -1)  # use cat instead of zcat if not gzipped

if [[ "$first_header" =~ ^@.+:[0-9]+:[A-Za-z0-9]+:[0-9]+:[0-9]+:[0-9]+:[0-9]+ ]]; then
    echo "CASAVA format detected"
    fastqc "$INPUT_FILE" -o "$OUTPUT_DIR" -t 6 --casava 
else
    echo "Not CASAVA format"
    fastqc "$INPUT_FILE" -o "$OUTPUT_DIR" -t 6 
fi



#frame the adapter args 
ADAPTER_ARGS= ""
for adapter_name in $ADAPTER_LIST; do
    adapter_seq="${!adapter_name}"
    ADAPTER_ARGS+=" -a $adapter_seq"
done



#running cutadapt to to remove adapters
echo "Removing adapters with cutadapt..."
if [[ ! -f "$INPUT_FILE" ]]; then
    echo "ERROR: Input FASTQ file $INPUT_FILE not found! Exiting..."
    exit 1
fi
echo "Current adapter: $ADAPTER_LIST"
cutadapt $ADAPTER_ARGS -m "$MINLEN" -M "$MAXLEN" -q "$MIN_QUALITY" -n 1 -O 5 --match-read-wildcards \
    -o trimmed_reads.fastq input_reads.fastq

echo "Cutadapt finished running!"

# Run FastQC on trimmed file
echo "Running FastQC on trimmed file..."
fastqc "$TRIMMED_FILE" -o "$FASTQC_OUTPUT_DIR" -t 6

echo "Running MINTmap..."
cd $MAINWORKDIR/MINT
./MINTmap.pl -f "$TRIMMED_FILE" -p "$MINT_OUTPUT_DIR"

# Report runtime
duration=$SECONDS
echo "Pipeline completed in $(($duration/60)) minutes and $(($duration % 60)) seconds."
echo "Results stored in:"
echo "  FastQC: $FASTQC_OUTPUT_DIR"
echo "  MINTmap: $MINT_OUTPUT_DIR"

