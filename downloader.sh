#!/bin/bash
set -e

# Configuration
WORK_DIR="/mnt/d/bioinformatics/SRA/tmp"
SRA_DIR="$WORK_DIR/sra"              # Directory to store downloaded SRA files
ACCESSION_FILE="$WORK_DIR/GSE156124.txt"  # File containing SRA accession numbers (one per line)
FASTQ_DIR="$WORK_DIR/fastq" 
# Create directories
mkdir -p "$SRA_DIR"

# Log file
LOG_FILE="$WORK_DIR/download_log1.txt"
echo "Starting SRA download and conversion at $(date)" > "$LOG_FILE"

# Check if accession file exists
if [[ ! -f "$ACCESSION_FILE" ]]; then
    echo "ERROR: Accession file $ACCESSION_FILE not found! Exiting..." | tee -a "$LOG_FILE"
    exit 1
fi

# Clean the accession file to remove trailing spaces
TMP_ACCESSION_FILE="$WORK_DIR/clean_GSE156124.txt"
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
    echo "Completed download and conversion for $ACCESSION" | tee -a "$LOG_FILE"

done < "$TMP_ACCESSION_FILE"
# Clean up temporary file
rm -f "$TMP_ACCESSION_FILE"

echo "all accessions downloaded"
EOF