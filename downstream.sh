#!/bin/bash

set -euo pipefail
IFS=$'\n\t'

timestamp() {
  date +"%Y-%m-%d %H:%M:%S"
}

# Accept arguments: metadata file (-r), mint folder path (-m)
RUNS_FILE=""
MINT_FOLDER_PATH=""

while getopts ":r:m:" opt; do
  case $opt in
    r) RUNS_FILE="$OPTARG" ;;
    m) MINT_FOLDER_PATH="$OPTARG" ;;
    \?) echo "Invalid option -$OPTARG" >&2; exit 1 ;;
    :) echo "Option -$OPTARG requires an argument." >&2; exit 1 ;;
  esac
done

if [[ -z "$RUNS_FILE" || -z "$MINT_FOLDER_PATH" ]]; then
  echo "Usage: $0 -r <runs.csv> -m <path/to/mint_folder>" >&2
  exit 1
fi

echo "[$(timestamp)] Starting script"
echo "Metadata file: $RUNS_FILE"
echo "MINT files folder: $MINT_FOLDER_PATH"

# Extract folder name from MINT_FOLDER_PATH
MINT_FOLDER_NAME=$(basename "$MINT_FOLDER_PATH")

echo "[$(timestamp)] Running MINTfilesmerger.py to merge MINT files..."
# Run the merger script and capture the merged file path from its stdout
MERGED_MINT_FILE=$(python3 MINTfilesmerger.py "$MINT_FOLDER_NAME" | grep "Merged file saved to:" | awk '{print $NF}')

if [[ ! -f "$MERGED_MINT_FILE" ]]; then
  echo "Error: Merged MINT file not found at $MERGED_MINT_FILE" >&2
  exit 1
fi
echo "[$(timestamp)] MINT files merged into: $MERGED_MINT_FILE"

echo "[$(timestamp)] Running meta_to_col.py..."
COL_DATA_DIR=$(python3 metaToCol.py "$RUNS_FILE")
if [[ $? -ne 0 ]]; then
  echo "Error: meta_to_col.py failed." >&2
  exit 1
fi
echo "[$(timestamp)] meta_to_col.py completed. colData saved in: $COL_DATA_DIR"

COL_DATA_FILE="${COL_DATA_DIR}/$(basename "${RUNS_FILE%.*}")_colData.csv"
if [[ ! -f "$COL_DATA_FILE" ]]; then
  echo "Error: Expected colData file not found at $COL_DATA_FILE" >&2
  exit 1
fi

echo "[$(timestamp)] Running diff_exp.R..."

SCRIPT_DIR="/raid/anirudh/bioinformatics/RNAseq_pipeline/script/trna_search"

diffexpout=$(docker run --rm \
  -v "$SCRIPT_DIR":/scripts \
  -v "$(dirname "$MERGED_MINT_FILE")":/data1 \
  -v "$(dirname "$COL_DATA_FILE")":/data2 \
  pegi3s/r_deseq2 \
  Rscript /scripts/diff_exp.R /data1/$(basename "$MERGED_MINT_FILE") /data2/$(basename "$COL_DATA_FILE") 2>&1)

if [[ $? -ne 0 ]]; then
  echo "$diffexpout"
  echo "Error: diff_exp.R failed." >&2
  exit 1
fi

# Extract file paths after marker line
file_list=$(echo "$diffexpout" | awk '/^### DESeq2 output files ###$/ {flag=1; next} flag && /^\/.*/')

echo "DESeq2 result files:"
echo "$file_list"


echo "[$(timestamp)] diff_exp.R completed successfully."
echo "[$(timestamp)] Script finished."


while IFS= read -r filepath; do
  # Trim leading/trailing whitespace from filepath
  filepath_trimmed=$(echo "$filepath" | xargs)
  if [[ -z "$filepath_trimmed" ]]; then
    continue
  fi
  echo "Filtering and merging DESeq2 results for $filepath_trimmed"
  python3 filter_deseq_results_and_merge.py "$MERGED_MINT_FILE" "$filepath_trimmed"
done <<< "$file_list"
