#!/bin/bash

# Check for required arguments
if [ $# -lt 2 ]; then
    echo "Usage: $0 <parameters_file> <accession_file>"
    echo "Example: $0 parameters.txt run_numbers.txt"
    exit 1
fi

PARAMETERS_FILE=$1
ACCESSION_FILE=$2
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


WORK_DIR="$MAINWORKDIR"
SRA_DIR="$WORK_DIR/sra"
FASTQ_DIR="$WORK_DIR/fastq"


mkdir -p "$SRA_DIR" "$FASTQ_DIR"

LOG_FILE="$WORK_DIR/download_log.txt"
echo "Starting SRA download and conversion at $(date)" > "$LOG_FILE"

TMP_ACCESSION_FILE="$WORK_DIR/RNAseq_pipeline/scripts/trna-search/clean_accessions.txt"
sed 's/[[:space:]]*$//' "$ACCESSION_FILE" > "$TMP_ACCESSION_FILE"

ORIG_DIR=$(pwd)

while IFS= read -r ACCESSION || [[ -n "$ACCESSION" ]]; do
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
    if ! fasterq-dump "$SRA_DIR/$ACCESSION.sra" \
        --outdir "$FASTQ_DIR/$ACCESSION" \
        --temp "$FASTQ_DIR/$ACCESSION/tmp" \
        --threads 4; then
        echo "ERROR: Failed to convert $ACCESSION, continuing..." | tee -a "$LOG_FILE"
        continue
    fi

    echo "Completed download and conversion for $ACCESSION" | tee -a "$LOG_FILE"

    echo "Starting tsRNA analysis pipeline for $ACCESSION..." | tee -a "$LOG_FILE"

    cd $MAINWORKDIR/RNAseq_pipeline/scripts/trna-search


    FASTQ_FILES=("$FASTQ_DIR/$ACCESSION"/*.fastq)
    if [[ ${#FASTQ_FILES[@]} -eq 0 ]]; then
        echo "ERROR: No FASTQ files found for $ACCESSION" | tee -a "$LOG_FILE"
        cd "$ORIG_DIR"
        continue
    fi

    if [[ ${#FASTQ_FILES[@]} -eq 1 ]]; then
        ./tRNA_pipeline.sh "${FASTQ_FILES[0]}" "$PARAMETERS_FILE" "$ACCESSION" || \
            echo "ERROR: Pipeline failed for $ACCESSION" | tee -a "$LOG_FILE"
    elif [[ ${#FASTQ_FILES[@]} -ge 2 ]]; then
        echo "Found paired-end FASTQ files: ${FASTQ_FILES[0]}, ${FASTQ_FILES[1]}" | tee -a "$LOG_FILE"
        ./Paired_tRNA_pipeline.sh "${FASTQ_FILES[0]}" "${FASTQ_FILES[1]}" "$PARAMETERS_FILE" "$ACCESSION" || \
            echo "ERROR: Pipeline failed for $ACCESSION" | tee -a "$LOG_FILE"
    else
        echo "ERROR: Unexpected number of FASTQ files for $ACCESSION" | tee -a "$LOG_FILE"
    fi

    cd "$ORIG_DIR"
done < "$TMP_ACCESSION_FILE"

echo "All SRA files processed." | tee -a "$LOG_FILE"
rm -f "$TMP_ACCESSION_FILE"
