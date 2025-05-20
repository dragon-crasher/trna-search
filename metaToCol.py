import pandas as pd
import sys
import os

def metaToCol(df, output_path):
    """
    Convert metadata to colData format and save to CSV.

    Parameters:
    - df: pandas DataFrame containing the metadata.
    - output_path: Path to save the filtered colData CSV.
    """
    # Select relevant columns and make a copy to avoid SettingWithCopyWarning
    filtered_df = df.loc[:, ['SRR_ID', 'Treatment']].copy()

    # Rename columns
    filtered_df.rename(columns={'SRR_ID': 'Sample', 'Treatment': 'Condition'}, inplace=True)

    # Save to CSV with Sample as the index (common for colData)
    filtered_df.set_index('Sample', inplace=True)
    filtered_df.sort_values(by='Condition', inplace=True)
    filtered_df.to_csv(output_path)

    print(f"Filtered metadata written to {output_path}", file=sys.stderr)  # Info to stderr

def main():
    if len(sys.argv) != 2:
        print("Usage: python metaToCol.py <metadata.csv>", file=sys.stderr)
        sys.exit(1)

    metadata_file = sys.argv[1]

    # Extract base name without extension
    name = os.path.splitext(os.path.basename(metadata_file))[0]

    try:
        df = pd.read_csv(metadata_file)
    except FileNotFoundError:
        print(f"Error: File '{metadata_file}' not found.", file=sys.stderr)
        sys.exit(1)
    except pd.errors.EmptyDataError:
        print(f"Error: File '{metadata_file}' is empty or invalid.", file=sys.stderr)
        sys.exit(1)

    # Check if required columns exist
    required_cols = {'SRR_ID', 'Treatment'}
    if not required_cols.issubset(df.columns):
        print(f"Error: Input file must contain columns: {required_cols}", file=sys.stderr)
        sys.exit(1)

    # Define output directory and path
    output_dir = "/raid/anirudh/bioinformatics/RNAseq_pipeline/diffexp"
    colDataPath = os.path.join(output_dir, f"{name}_colData.csv")

    # Run conversion and save
    metaToCol(df, colDataPath)

    # Print output directory path to stdout for bash script capture
    print(output_dir)

if __name__ == "__main__":
    main()
