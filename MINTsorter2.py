import os
import re
import shutil

def process_files(folder_path):
    # Step 1: Rename files containing 'counts'
    for filename in os.listdir(folder_path):
        if filename.endswith('.txt') and 'expression' in filename:
            accession_code = re.split(r"[_-]", filename)[0]
            if 'exclusive' in filename:
                suffix = 'e'
            elif 'ambiguous' in filename:
                suffix = 'a'
            else:
                continue
            selected_folder = os.path.join(folder_path, "selected")
            os.makedirs(selected_folder, exist_ok=True)
            new_filename = f"{suffix}_{accession_code}.txt"
            src_path = os.path.join(folder_path, filename)
            dst_path = os.path.join(selected_folder, new_filename)
            print(f"Moving: {src_path} -> {dst_path}")
            shutil.move(src_path, dst_path)
    
    # Step 2: Move all files containing 'selected' to the 'selected' folder
    
               

# Example usage:
process_files("/mnt/d/bioinformatics/MINT/outputs/")
