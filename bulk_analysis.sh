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
  echo "Usage: $0 -p <parameters.txt> -d <Run_numbers.txt> -a <run cutadapt? :1|0>" >&2
  exit 1
fi

echo "Parameters file: $PARAMETERS_FILE"
echo "Accession File: $ACCESSION_FILE"
echo "Run Cutadapt: $CUTADAPT_RUN"

source "$PARAMETERS_FILE"


# Define WORK_DIR before LOG_FILE
LOG_DIR="$MAINWORKDIR/RNAseq_pipeline/script/trna_search/logs"
mkdir -p "$LOG_DIR"

# Set log filename based on accession file basename (without extension)
ACCESSION_BASENAME=$(basename "$ACCESSION_FILE")
LOG_BASENAME="${ACCESSION_BASENAME%.*}"
LOG_FILE="$LOG_DIR/${LOG_BASENAME}_download_log.txt"

# Initialize log file
echo "$(date '+%Y-%m-%d %H:%M:%S') - Starting SRA download and conversion" > "$LOG_FILE"

log() {
    echo "$(date '+%Y-%m-%d %H:%M:%S') - $*" | tee -a "$LOG_FILE"
}

# Validate input files
if [[ ! -f "$PARAMETERS_FILE" ]]; then
    log "ERROR: Parameters file $PARAMETERS_FILE not found! Exiting..."
    exit 1
fi

if [[ ! -f "$ACCESSION_FILE" ]]; then
    log "ERROR: Accession file $ACCESSION_FILE not found! Exiting..."
    exit 1
fi

WORK_DIR="$MAINWORKDIR/SRA"
SRA_DIR="$WORK_DIR/sra"
FASTQ_DIR="$WORK_DIR/fastq"

mkdir -p "$SRA_DIR" "$FASTQ_DIR"

TMP_ACCESSION_FILE="$MAINWORKDIR/RNAseq_pipeline/script/trna_search/clean_accessions.txt"
sed 's/[[:space:]]*$//' "$ACCESSION_FILE" > "$TMP_ACCESSION_FILE"

TOTAL_ACCESSIONS=$(grep -v -e '^\s*$' -e '^\s*#' "$TMP_ACCESSION_FILE" | wc -l)
COUNT=0
BAR_WIDTH=40

print_progress_bar() {
  local progress=$1
  local total=$2
  local width=$3

  local percent=$(( progress * 100 / total ))
  local filled=$(( progress * width / total ))
  local empty=$(( width - filled ))

  # Construct bar string
  local bar_filled=$(printf "%${filled}s" | tr ' ' '#')
  local bar_empty=$(printf "%${empty}s")

  # Print progress bar with percentage
  echo -ne "Progress: [${bar_filled}${bar_empty}] ${percent}% ($progress/$total) \r"
}



ORIG_DIR=$(pwd)

# Define number of CPU threads (adjust as needed)
cputhreads=4

while IFS= read -r ACCESSION || [[ -n "$ACCESSION" ]]; do
    
    SECONDS=0
    
    [[ -z "$ACCESSION" || "$ACCESSION" == \#* ]] && continue
    ACCESSION=$(echo "$ACCESSION" | xargs)

    log "Processing accession: $ACCESSION"
    mkdir -p "$FASTQ_DIR/$ACCESSION"

    log "Downloading $ACCESSION..."
    if ! prefetch "$ACCESSION" --output-directory "$SRA_DIR" >> "$LOG_FILE" 2>&1; then
        log "ERROR: Failed to download $ACCESSION, continuing..."
        continue
    fi

    log "Converting $ACCESSION to FASTQ..."
    if ! fasterq-dump "$SRA_DIR/$ACCESSION/$ACCESSION.sra" \
        --outdir "$FASTQ_DIR/$ACCESSION" \
        --temp "$FASTQ_DIR/$ACCESSION/tmp" \
        --threads "$cputhreads" \
        --force >> "$LOG_FILE" 2>&1; then
        log "ERROR: Failed to convert $ACCESSION, continuing..."
    fi

    log "Completed download and conversion for $ACCESSION"

    log "Starting tsRNA analysis pipeline for $ACCESSION..."

    cd "$MAINWORKDIR/RNAseq_pipeline/script/trna_search" || { log "ERROR: Failed to cd to script directory"; exit 1; }

    FASTQ_FILES=("$FASTQ_DIR/$ACCESSION"/*.fastq)
    if [[ ${#FASTQ_FILES[@]} -eq 0 ]]; then
        log "ERROR: No FASTQ files found for $ACCESSION"
        cd "$ORIG_DIR" || exit
        continue
    fi

    if [[ ${#FASTQ_FILES[@]} -eq 1 ]]; then
        ./tRNA_pipeline.sh "${FASTQ_FILES[0]}" "$PARAMETERS_FILE" "$ACCESSION" "$CUTADAPT_RUN" || \
        log "ERROR: Pipeline failed for $ACCESSION"
    elif [[ ${#FASTQ_FILES[@]} -ge 2 ]]; then
        log "Found paired-end FASTQ files: ${FASTQ_FILES[0]}, ${FASTQ_FILES[1]}"
        ./Paired_tRNA_pipeline.sh "${FASTQ_FILES[0]}" "${FASTQ_FILES[1]}" "$PARAMETERS_FILE" "$ACCESSION" "$CUTADAPT_RUN" || \
        log "ERROR: Pipeline failed for $ACCESSION"
    else
        log "ERROR: Unexpected number of FASTQ files for $ACCESSION"
    fi

    duration=$SECONDS
    log "Pipeline for $ACCESSION completed in $(($duration / 60)) minutes and $(($duration % 60)) seconds."

    COUNT=$((COUNT + 1))
    print_progress_bar "$COUNT" "$TOTAL_ACCESSIONS" "$BAR_WIDTH"

    cd "$ORIG_DIR" || exit
done < "$TMP_ACCESSION_FILE"

# python3 MINTsorter2.py "$PROJECT_NAME"
log "All SRA files processed."

rm -f "$TMP_ACCESSION_FILE"
