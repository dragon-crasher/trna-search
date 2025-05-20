import pandas as pd
import sys

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

    print(f"Filtered metadata written to {output_path}")

def main():
    if len(sys.argv) != 3:
        print("Usage: python metaToCol.py <metadata.csv> <output_metadata.csv>")
        sys.exit(1)

    metadata_file = sys.argv[1]
    output_metadata_file = sys.argv[2]

    try:
        df = pd.read_csv(metadata_file)
    except FileNotFoundError:
        print(f"Error: File '{metadata_file}' not found.")
        sys.exit(1)
    except pd.errors.EmptyDataError:
        print(f"Error: File '{metadata_file}' is empty or invalid.")
        sys.exit(1)

    # Check if required columns exist
    required_cols = {'SRR_ID', 'Treatment'}
    if not required_cols.issubset(df.columns):
        print(f"Error: Input file must contain columns: {required_cols}")
        sys.exit(1)

    metaToCol(df, output_metadata_file)

if __name__ == "__main__":
    main()
