#!/bin/bash
set -e


# Check for required arguments
if [ $# -lt 2 ]; then
    echo "Usage: $0 <parameters_file> <accession_file>"
    echo "Example: $0 parameters.txt run_numbers.txt"
    exit 1
fi

PARAMETERS_FILE=$1
ACCESSION_FILE=${2:-parameters.txt}

# Validate input files
if [[ ! -f "$PARAMETERS_FILE" ]]; then
    echo "ERROR: Parameters file $PARAMETERS_FILE not found! Exiting..."
    exit 1
fi

if [[ ! -f "$ACCESSION_FILE" ]]; then
    echo "ERROR: Accession file $ACCESSION_FILE not found! Exiting..."
    exit 1
fi

# Configuration
WORK_DIR="/mnt/d/bioinformatics/SRA/tmp"
SRA_DIR="$WORK_DIR/sra"              # Directory to store downloaded SRA files
FASTQ_DIR="$WORK_DIR/fastq"          # Directory to store converted FASTQ files

# Create directories
mkdir -p "$SRA_DIR" "$FASTQ_DIR"

# Log file
LOG_FILE="$WORK_DIR/download_log.txt"
echo "Starting SRA download and conversion at $(date)" > "$LOG_FILE"

# Clean the accession file to remove trailing spaces
TMP_ACCESSION_FILE="$WORK_DIR/clean_accessions.txt"
sed 's/[[:space:]]*$//' "$ACCESSION_FILE" > "$TMP_ACCESSION_FILE"   

# Process all accessions
while IFS= read -r ACCESSION || [[ -n "$ACCESSION" ]]; do
    # Skip empty lines and comments
    [[ -z "$ACCESSION" || "$ACCESSION" == \#* ]] && continue
    
    # Trim any leading or trailing whitespace
    ACCESSION=$(echo "$ACCESSION" | xargs)
    
    echo "Processing accession: $ACCESSION" | tee -a "$LOG_FILE"
    
    # Create output directory for this accession
    mkdir -p "$FASTQ_DIR/$ACCESSION"
    
    # Step 1: Download SRA file using prefetch
    echo "Downloading $ACCESSION..." | tee -a "$LOG_FILE"
    prefetch "$ACCESSION" --output-directory "$SRA_DIR" || {
        echo "ERROR: Failed to download $ACCESSION, continuing with next accession..." | tee -a "$LOG_FILE"
        continue
    }
    
    # Step 2: Convert SRA to FASTQ using fasterq-dump
    echo "Converting $ACCESSION to FASTQ..." | tee -a "$LOG_FILE"
    fasterq-dump "$SRA_DIR/$ACCESSION/$ACCESSION.sra" \
        --outdir "$FASTQ_DIR/$ACCESSION" \
        --temp "$FASTQ_DIR/$ACCESSION/tmp" \
        --threads 4 || {
        echo "ERROR: Failed to convert $ACCESSION, continuing with next accession..." | tee -a "$LOG_FILE"
        continue
    }
    
    echo "Completed download and conversion for $ACCESSION" | tee -a "$LOG_FILE"
    
    # Run tsRNA analysis pipeline for this accession
    echo "Starting tsRNA analysis pipeline for $ACCESSION..." | tee -a "$LOG_FILE"
    
    # Change directory to where the pipeline script is located
    cd /mnt/d/bioinformatics/RNAseq_pipeline/script
    
    # Find the FASTQ file(s) for this accession
    FASTQ_FILES=("$FASTQ_DIR/$ACCESSION"/*.fastq)
    
    if [[ ${#FASTQ_FILES[@]} -eq 0 ]]; then
        echo "ERROR: No FASTQ files found for $ACCESSION" | tee -a "$LOG_FILE"
        continue
    fi
    
    # If single-end data (one FASTQ file)
    if [[ ${#FASTQ_FILES[@]} -eq 1 ]]; then
        # Run the pipeline script for this accession
        ./tRNA_pipeline.sh "${FASTQ_FILES[0]}" "$PARAMETERS_FILE" "$ACCESSION_FILE" "$ACCESSION" || {
            echo "ERROR: Pipeline failed for $ACCESSION" | tee -a "$LOG_FILE"
        }
    else
        # If paired-end data, you may need to modify your pipeline script or handle differently
        echo "Found multiple FASTQ files for $ACCESSION. Using first file: ${FASTQ_FILES[0]}" | tee -a "$LOG_FILE"
        ./tRNA_pipeline.sh "${FASTQ_FILES[0]}" "$PARAMETERS_FILE" "$ACCESSION_FILE" "$ACCESSION" || {
            echo "ERROR: Pipeline failed for $ACCESSION" | tee -a "$LOG_FILE"
        }
    fi
    
    # Go back to the original directory for the next iteration
    cd - > /dev/null
    
done < "$TMP_ACCESSION_FILE"

echo "All SRA files have been downloaded, converted to FASTQ, and analyzed" | tee -a "$LOG_FILE"

# Clean up temporary file
rm -f "$TMP_ACCESSION_FILE"