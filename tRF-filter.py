import pandas as pd
import os
import sys

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

def parse_parameters_file(filepath):
    params = {}
    with open(filepath, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue  # skip empty lines and comments
            if '=' in line:
                key, value = line.split('=', 1)
                key = key.strip()
                value = value.strip().strip('"').strip("'")  # remove quotes if any
                # Convert numeric values if possible
                if value.isdigit():
                    value = int(value)
                params[key] = value
    return params

if __name__ == "__main__":
    # Check if the script is run with the correct number of arguments

    if len(sys.argv) < 2:
        print("Usage: python tRF-filter.py [parameters_file.txt] [tRF_fragments_folder_path] [sample_names]")
        sys.exit(1)

    # import the parameters file and folder path from command line arguments
    parameters_file = sys.argv[1]

    # import the parameters to dictoionary format
    parameters = parse_parameters_file(parameters_file)

    folder_path = sys.argv[2]

    if not os.path.exists(folder_path):
        print(f"Error: The folder path '{folder_path}' does not exist.")
        sys.exit(1)
    if not os.path.isdir(folder_path):
        print(f"Error: The path '{folder_path}' is not a directory.")
        sys.exit(1)
   


    # define the output folder path
    out_folder_path = parameters.get("output_path", "/mnt/d/bioinformatics/RNAseq_pipeline/diffexp/")

    dir_path = out_folder_path
    file_count = len([entry for entry in os.listdir(dir_path) if os.path.isfile(os.path.join(dir_path, entry))]) + 1

    base_name = sys.argv[3] if len(sys.argv) > 3 else str(file_count)

    csv_output_name = os.path.join(out_folder_path, f"{base_name}_merged.csv")
    trf_sequences_name = os.path.join(out_folder_path, f"{base_name}_trf_sequences.csv")
    
           
    read_and_merge_dataframes(
        folder_path,
        count_threshold=0,
        output_csv= csv_output_name,
        sequence_output_csv= trf_sequences_name
    )
