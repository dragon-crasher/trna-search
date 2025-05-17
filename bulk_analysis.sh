#!/bin/bash

PARAMETERS_FILE=""
ACCESSION_FILE=""
CUTADAPT_RUN=1

while getopts ":p:d:a:" opt; do
  case $opt in
    p) PARAMETERS_FILE="$OPTARG" ;;
    d) ACCESSION_FILE="$OPTARG" ;;
    a) CUTADAPT_RUN="$OPTARG" ;;
    \?) echo "Invalid option -$OPTARG" >&2; exit 1 ;;
    :) echo "Option -$OPTARG requires an argument." >&2; exit 1 ;;
  esac
done

# Check for required arguments
if [[ -z "$PARAMETERS_FILE" || -z "$ACCESSION_FILE" || -z "$CUTADAPT_RUN" ]]; then
  echo "Usage: $0 -p <parameters.txt> -d <Run_numbers.txt> -a <run cutadapt? :1|0>
  "
  exit 1
fi

echo "Parameters file: $PARAMETERS_FILE"
echo "Accession File: $ACCESSION_FILE"
echo "Run Cutadapt: $CUTADAPT_RUN"


source "$PARAMETERS_FILE"
# Validate input files
if [[ ! -f "$PARAMETERS_FILE" ]]; then
    echo "ERROR: Parameters file $PARAMETERS_FILE not found! Exiting..."
    exit 1
fi

if [[ ! -f "$ACCESSION_FILE" ]]; then
    echo "ERROR: Accession file $ACCESSION_FILE not found! Exiting..."
    exit 1
fi


WORK_DIR="$MAINWORKDIR/SRA"
SRA_DIR="$WORK_DIR/sra"
FASTQ_DIR="$WORK_DIR/fastq"


mkdir -p "$SRA_DIR" "$FASTQ_DIR"

LOG_FILE="$WORK_DIR/download_log.txt"
echo "Starting SRA download and conversion at $(date)" > "$LOG_FILE"

TMP_ACCESSION_FILE="$MAINWORKDIR/RNAseq_pipeline/script/trna_search/clean_accessions.txt"
sed 's/[[:space:]]*$//' "$ACCESSION_FILE" > "$TMP_ACCESSION_FILE"

ORIG_DIR=$(pwd)

while IFS= read -r ACCESSION || [[ -n "$ACCESSION" ]]; do
    
    SECONDS=0

    [[ -z "$ACCESSION" || "$ACCESSION" == \#* ]] && continue
    ACCESSION=$(echo "$ACCESSION" | xargs)

    echo "Processing accession: $ACCESSION" | tee -a "$LOG_FILE"
    mkdir -p "$FASTQ_DIR/$ACCESSION"

    echo "Downloading $ACCESSION..." | tee -a "$LOG_FILE"
    if ! prefetch "$ACCESSION" --output-directory "$SRA_DIR"; then
        echo "ERROR: Failed to download $ACCESSION, continuing..." | tee -a "$LOG_FILE"
        continue
    fi

    echo "Converting $ACCESSION to FASTQ..." | tee -a "$LOG_FILE"
    if ! fasterq-dump "$SRA_DIR/$ACCESSION/$ACCESSION.sra" \
        --outdir "$FASTQ_DIR/$ACCESSION" \
        --temp "$FASTQ_DIR/$ACCESSION/tmp" \
        --threads 12 \
        --force; then
        echo "ERROR: Failed to convert $ACCESSION, continuing..." | tee -a "$LOG_FILE"
    fi

    echo "Completed download and conversion for $ACCESSION" | tee -a "$LOG_FILE"

    echo "Starting tsRNA analysis pipeline for $ACCESSION..." | tee -a "$LOG_FILE"

    cd $MAINWORKDIR/RNAseq_pipeline/script/trna_search


    FASTQ_FILES=("$FASTQ_DIR/$ACCESSION"/*.fastq)
    if [[ ${#FASTQ_FILES[@]} -eq 0 ]]; then
        echo "ERROR: No FASTQ files found for $ACCESSION" | tee -a "$LOG_FILE"
        cd "$ORIG_DIR"
        continue
    fi

    if [[ ${#FASTQ_FILES[@]} -eq 1 ]]; then
        ./tRNA_pipeline.sh "${FASTQ_FILES[0]}" "$PARAMETERS_FILE" "$ACCESSION" "$CUTADAPT_RUN" || \
        echo "ERROR: Pipeline failed for $ACCESSION" | tee -a "$LOG_FILE"
    elif [[ ${#FASTQ_FILES[@]} -ge 2 ]]; then
        echo "Found paired-end FASTQ files: ${FASTQ_FILES[0]}, ${FASTQ_FILES[1]}" | tee -a "$LOG_FILE"
        ./Paired_tRNA_pipeline.sh "${FASTQ_FILES[0]}" "${FASTQ_FILES[1]}" "$PARAMETERS_FILE" "$ACCESSION" "$CUTADAPT_RUN" || \
        echo "ERROR: Pipeline failed for $ACCESSION" | tee -a "$LOG_FILE"
    else
        echo "ERROR: Unexpected number of FASTQ files for $ACCESSION" | tee -a "$LOG_FILE"
    fi


    duration=$SECONDS
    echo "Pipeline for $ACCESSION completed in $(($duration / 60)) minutes and $(($duration % 60)) seconds."
    
    cd "$ORIG_DIR"
done < "$TMP_ACCESSION_FILE"

echo "All SRA files processed." | tee -a "$LOG_FILE"
rm -f "$TMP_ACCESSION_FILE"
