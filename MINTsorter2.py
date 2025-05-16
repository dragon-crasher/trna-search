import os
import re
import shutil

def process_files(folder_path):
    # Step 1: Rename files containing 'counts'
    for filename in os.listdir(folder_path):
        if filename.endswith('.txt') and 'counts' in filename:
            accession_code = re.split(r"[_-]", filename)[0]
            if 'exclusive' in filename:
                suffix = 'exclusive'
            elif 'ambiguous' in filename:
                suffix = 'ambiguous'
            else:
                continue
            new_filename = f"{accession_code}_selected_{suffix}.txt"
            src_path = os.path.join(folder_path, filename)
            dst_path = os.path.join(folder_path, new_filename)
            print(f"Renaming: {src_path} -> {dst_path}")
            os.rename(src_path, dst_path)
    
    # Step 2: Move all files containing 'selected' to the 'selected' folder
    selected_folder = os.path.join(folder_path, "selected")
    os.makedirs(selected_folder, exist_ok=True)
    for filename in os.listdir(folder_path):
        if filename.endswith('.txt') and 'selected' in filename:
            src_path = os.path.join(folder_path, filename)
            dst_path = os.path.join(selected_folder, filename)
            print(f"Moving: {src_path} -> {dst_path}")
            shutil.move(src_path, dst_path)

# Example usage:
process_files("/mnt/d/bioinformatics/MINT/outputs/")
