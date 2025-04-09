import pandas as pd
import os

def read_and_merge_dataframes(folder_path, count_threshold=10, output_csv="merged_data.csv", sequence_output_csv="trf_sequences.csv"):
    merged_df = None
    sample_names = []
    trf_sequences_df = None
    
    for filename in os.listdir(folder_path):
        file_path = os.path.join(folder_path, filename)
        try:
            df = pd.read_csv(file_path, sep='\t')
            filename_without_extension = os.path.splitext(filename)[0]

            required_columns = ['tRF sequence', 'Unnormalized read counts', 'MINTbase Unique ID']  
            if all(col in df.columns for col in required_columns):
                # Create a separate DataFrame for tRF sequences and IDs
                trf_sequences = df[['MINTbase Unique ID', 'tRF sequence']].copy()
                
                if trf_sequences_df is None:
                    trf_sequences_df = trf_sequences
                else:
                    trf_sequences_df = pd.concat([trf_sequences_df, trf_sequences], ignore_index=True)

                # Filter and merge data
                df = df[['tRF sequence', 'Unnormalized read counts']].copy()
                df = df[df['Unnormalized read counts'] >= count_threshold]
                df.rename(columns={'Unnormalized read counts': filename_without_extension}, inplace=True)

                if merged_df is None:
                    merged_df = df
                else:
                    merged_df = pd.merge(merged_df, df, on='tRF sequence', how='outer')

                sample_names.append(filename_without_extension)

            else:
                print(f"Warning: Required columns not found in {filename}. Skipping.")
                continue

        except Exception as e:
            print(f"Error reading {filename}: {e}")

    if merged_df is not None:
        merged_df.fillna(0, inplace=True)
        merged_df.set_index('tRF sequence', inplace=True)  
        merged_df.to_csv(output_csv)  
        
        # Export tRF sequences with their IDs
        if trf_sequences_df is not None:
            trf_sequences_df.drop_duplicates(inplace=True)
            trf_sequences_df.to_csv(sequence_output_csv, index=False)
            print(f"tRF sequences exported to {sequence_output_csv}")
        
        print(f"Merged DataFrame exported to {output_csv}")
    else:
        print("No DataFrames were merged.")

if __name__ == "__main__":
    folder_path = "/mnt/d/bioinformatics/MINT/MINTmap-release-v1.0/MINTmap-release-v1.0/outputs/CLL/chosen"  
    read_and_merge_dataframes(
        folder_path,
        count_threshold=0,
        output_csv="/mnt/d/bioinformatics/RNAseq_pipeline/diffexp/merged_tRF_data.csv",
        sequence_output_csv="/mnt/d/bioinformatics/RNAseq_pipeline/diffexp/trf_sequences.csv"
    )
