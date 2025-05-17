import os
import shutil
import re

def rename_and_move_counts_files(parent_dir):
    # Walk through each subdirectory in the parent directory
    for root, dirs, files in os.walk(parent_dir):
        # Skip the parent directory itself
        if root == parent_dir:
            continue
        folder_name = os.path.basename(root)
        for file in files:
            if "expression" in file:
                old_file_path = os.path.join(root, file)
                # Determine the new name based on file content
                if "exclusive" in file:
                    new_file_name = f"e_{folder_name}"
                elif "ambiguous" in file:
                    new_file_name = f"a_{folder_name}"
                else:
                    continue  # Skip files that don't match either
                # Preserve the file extension
                _, ext = os.path.splitext(file)
                new_file_name += ext
                new_file_path = os.path.join(root, new_file_name)
                # Rename the file
                os.rename(old_file_path, new_file_path)
                print(f"Renamed: {old_file_path} -> {new_file_path}")

                # Prepare the 'selected' folder path
                selected_folder = os.path.join(root, "selected")
                os.makedirs(selected_folder, exist_ok=True)
                # Move the renamed file to the 'selected' folder
                final_path = os.path.join(selected_folder, new_file_name)
                shutil.move(new_file_path, final_path)
                print(f"Moved: {new_file_path} -> {final_path}")

# Example usage:
rename_and_move_counts_files("/mnt/d/bioinformatics/MINT/outputs/")


def move_selected_files_to_results(parent_dir, results_dir):
    # Ensure the results directory exists
    os.makedirs(results_dir, exist_ok=True)

    # Iterate over all items in the parent directory
    for item in os.listdir(parent_dir):
        item_path = os.path.join(parent_dir, item)
        # Check if the item is a directory (accession folder)
        if os.path.isdir(item_path):
            selected_path = os.path.join(item_path, "selected")
            # Check if the 'selected' subfolder exists
            if os.path.isdir(selected_path):
                for filename in os.listdir(selected_path):
                    if filename.endswith('.txt'):
                        src_file = os.path.join(selected_path, filename)
                        dst_file = os.path.join(results_dir, filename)
                        # Only move files (not subdirectories)
                        if os.path.isfile(src_file):
                            print(f"Moving {src_file} -> {dst_file}")
                            shutil.move(src_file, dst_file)

# Usage

move_selected_files_to_results(
    "/mnt/d/bioinformatics/MINT/outputs",
    "/mnt/d/bioinformatics/MINT/outputs/selected"
)



def standardize_filenames(folder_path):
    for filename in os.listdir(folder_path):
        if not filename.endswith('.txt'):
            continue

        # Try to find the accession number (SRR followed by digits)
        match = re.search(r'(SRR\d+)', filename)
        if not match:
            continue  # Skip files without an accession number

        accession = match.group(1)

        if 'exclusive' in filename:
            new_filename = f"e_{accession}.txt"
        elif 'ambiguous' in filename:
            new_filename = f"a_{accession}.txt"
        else:
            continue  # Skip files that are neither

        src = os.path.join(folder_path, filename)
        dst = os.path.join(folder_path, new_filename)
        print(f"Renaming: {filename} -> {new_filename}")
        os.rename(src, dst)

# Usage:
#rename_and_move_counts_files("/mnt/d/bioinformatics/MINT/outputs/")



