#!/bin/bash
set -euo pipefail

# Check for args and set defaults
if [ $# -lt 2 ]; then
    echo "Usage: $0 <input_fastq_file> <parameter_file> [output_dir] <run_cutadapt: 1|0>"
    echo "Example: $0 sample1.fastq parameters.txt sample1_results"
    exit 1
fi

INPUT_FILE=$1
PARAMETER_FILE=$2
run_cutadapt=${4:-1}

# Check if input files exist before sourcing
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

# Source parameters file
source "$PARAMETER_FILE"

# Set output directory, use third argument if provided, otherwise use basename of input file
if [ $# -ge 3 ]; then
    OUTPUT_DIR=$3
else
    OUTPUT_DIR=$(basename "$INPUT_FILE" .fastq)
fi

# Import list of adapters and parameters from adapters.txt
while IFS='=' read -r name seq; do
    # Skip empty lines or lines starting with #
    [[ -z "$name" || "$name" == \#* ]] && continue

    # Trim whitespace from name and seq
    name=$(echo "$name" | xargs)
    seq=$(echo "$seq" | xargs)

    # Remove surrounding quotes from seq if present
    seq=${seq#\"}
    seq=${seq%\"}

    # Declare variable dynamically
    declare "$name=$seq"
done < adapters.txt

# Source parameters again in case adapters.txt defines variables used in params
source "$PARAMETER_FILE"

echo "Processing file: $INPUT_FILE"
echo "Output directory: $OUTPUT_DIR"
echo "Current shell: $SHELL"


# Define output directories
MINT_OUTPUT_DIR="$MAINWORKDIR/MINT/outputs/$OUTPUT_DIR"
FASTQC_OUTPUT_DIR="$MAINWORKDIR/RNAseq_pipeline/data/fastqc/$OUTPUT_DIR/"

# Create output directories if they don't exist
mkdir -p "$MINT_OUTPUT_DIR"
mkdir -p "$FASTQC_OUTPUT_DIR"


# Define accession base name and trimmed file name
ACCESSION=$(basename "$INPUT_FILE" .fastq)
ACCESSION_trimmed="${ACCESSION}_trimmed.fastq"

# Perform quality check on input file
echo "Running FastQC on input file..."
first_header=$(head -1 "$INPUT_FILE")

if [[ "$first_header" =~ ^@.+:[0-9]+:[A-Za-z0-9]+:[0-9]+:[0-9]+:[0-9]+:[0-9]+ ]]; then
    echo "CASAVA format detected"
    fastqc "$INPUT_FILE" -o "$FASTQC_OUTPUT_DIR" -t $cputhreads --casava
else
    echo "Not CASAVA format"
    fastqc "$INPUT_FILE" -o "$FASTQC_OUTPUT_DIR" -t $cputhreads
fi

if ((run_cutadapt == 1)); then
    echo "Running cutadapt..."
    # Frame the adapter args
    ADAPTER_ARGS=""
    	for adapter_name in $ADAPTER3_LIST; do
        	if [[ -z "${!adapter_name:-}" ]]; then
        		echo "ERROR: Adapter sequence for '$adapter_name' is not set. Check adapters.txt."
        		exit 1
    		fi
		adapter_seq="${!adapter_name}"
        	ADAPTER_ARGS+=" -a $adapter_seq"
	done
	for adapt_name in $ADAPTER5_LIST; do
		if [[ -z "${!adapter_name:-}" ]]; then
        		echo "ERROR: Adapter sequence for '$adapter_name' is not set. Check adapters.txt."
        		exit 1
    		fi
		
		adapter_seq="${!adapter_name}"
		ADAPTER_ARGA+=" -g $adapter_seq"
    	done

    # Running cutadapt to remove adapters
    echo "Removing adapters with cutadapt..."
    cutadapt $ADAPTER_ARGS -m "$MINLEN" -M "$MAXLEN" -q "$MIN_QUALITY" -n 1 -j $cpucores -O 5 --match-read-wildcards \
        -o "$MAINWORKDIR/SRA/fastq/$ACCESSION/$ACCESSION_trimmed" "$INPUT_FILE"
    echo "Cutadapt finished running!"

    # Run FastQC on trimmed file
    echo "Running FastQC on trimmed file..."
    fastqc "$MAINWORKDIR/SRA/fastq/$ACCESSION/$ACCESSION_trimmed" -o "$FASTQC_OUTPUT_DIR" -t $cputhreads
    RUNNING=$MAINWORKDIR/SRA/fastq/$ACCESSION/$ACCESSION_trimmed
else
    echo "Skipping cutadapt step..."
    RUNNING=$INPUT_FILE
fi



# Run MINTmap
echo "Running MINTmap..."
cd "$MAINWORKDIR/MINT"

#version 1 of MINTmap
#./MINTmap.pl -f "$RUNNING" -p "$MINT_OUTPUT_DIR"

#version 2 of MINTmap
cd "$MAINWORKDIR/MINT/outputs/"
conda run -n mintmap38 MINTmap -p $OUTPUT_DIR "$RUNNING" 
# Report runtime

echo "Results stored in:"
echo "  FastQC: $FASTQC_OUTPUT_DIR"
echo "  MINTmap: $MINT_OUTPUT_DIR"
