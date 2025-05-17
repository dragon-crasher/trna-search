import os
import re
import sys
import shutil

def process_files(folder_path, name="selected"):
    selected_folder = os.path.join(folder_path, name)
    os.makedirs(selected_folder, exist_ok=True)

    for filename in os.listdir(folder_path):
        if filename.endswith('.txt'):  # Remove 'expression' filter if not needed
            accession_code = re.split(r"[_-]", filename)[0]
            if 'exclusive' in filename:
                suffix = 'e'
            elif 'ambiguous' in filename:
                suffix = 'a'
            else:
                continue

            # Check for R1 or R2 in the filename (case-insensitive)
            r1_match = re.search(r'R1', filename, re.IGNORECASE)
            r2_match = re.search(r'R2', filename, re.IGNORECASE)

            if r1_match:
                new_filename = f"{suffix}_{accession_code}_r1.txt"
            elif r2_match:
                new_filename = f"{suffix}_{accession_code}_r2.txt"
            else:
                new_filename = f"{suffix}_{accession_code}.txt"

            src_path = os.path.join(folder_path, filename)
            dst_path = os.path.join(selected_folder, new_filename)
            print(f"Moving: {src_path} -> {dst_path}")
            shutil.move(src_path, dst_path)

# Example usage:
if len(sys.argv) != 2:
    print("Usage: python MINTsorter2.py <name>")
    sys.exit(1)
name = sys.argv[1]
process_files("/mnt/d/bioinformatics/MINT/outputs/", name)
