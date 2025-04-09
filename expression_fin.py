import pandas as pd

def merge_csv_files(file1, file2, output_file):
    try:
        # Read the two CSV files
        df1 = pd.read_csv(file1)
        df2 = pd.read_csv(file2)
        
        # Check if 'tRF sequence' column exists in both files
        if 'tRF sequence' not in df1.columns or 'tRF sequence' not in df2.columns:
            raise ValueError("Both files must contain the 'tRF sequence' column.")
        
        # Merge the two DataFrames based on 'tRF sequence'
        merged_df = pd.merge(df1, df2, on='tRF sequence', how='outer')
        
        # Export the merged DataFrame to a new CSV file
        merged_df.to_csv(output_file, index=False)
        
        print(f"Merged file saved to {output_file}")
    except Exception as e:
        print(f"Error: {e}")

if __name__ == "__main__":
    # File paths
    file1 = "/path/to/first_file.csv"
    file2 = "/path/to/second_file.csv"
    output_file = "/path/to/merged_file.csv"
    
    # Merge the files
    merge_csv_files(file1, file2, output_file)
