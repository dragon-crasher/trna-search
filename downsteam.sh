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
COL_DATA_DIR=$(python3 meta_to_col.py "$RUNS_FILE")
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
Rscript diff_exp.R "$COL_DATA_FILE" "$MERGED_MINT_FILE"
if [[ $? -ne 0 ]]; then
  echo "Error: diff_exp.R failed." >&2
  exit 1
fi

echo "[$(timestamp)] diff_exp.R completed successfully."
echo "[$(timestamp)] Script finished."
