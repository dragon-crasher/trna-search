import pandas as pd
import sys
import os

if len(sys.argv) != 3:
    print("Usage: python filter_deseq_results_and_merge.py <merged_mint_files.csv> <deseq2_results.csv>")
    sys.exit(1)

# === USER INPUTS ===
other_file = sys.argv[1]
deseq2_results_file = sys.argv[2]

output_file = "/raid/anirudh/bioinformatics/RNAseq_pipeline/data/" + os.path.splitext(os.path.basename(deseq2_results_file))[0] + "_leads.csv"
significance_threshold = 0.05

# === 1. Read DESeq2 results ===
df = pd.read_csv(deseq2_results_file)
df.rename(columns={df.columns[0]: 'License Plate'}, inplace=True)

# === 2. Filter for significant results (padj not NA and < threshold) ===
df_filtered = df[df['padj'].notna() & (df['padj'] < significance_threshold)]

# === 3. Read the other file (first 3 columns only) ===
other_df = pd.read_csv(other_file)
other_df = other_df[['License Plate', 'tRF sequence', 'tRF type(s)']]

# === 4. Merge on 'MINTbase Unique ID' ===
merged = pd.merge(other_df, df_filtered, on='License Plate', how='inner')
merged = merged.sort_values(by = 'padj', ascending=True)
# === 5. Write output ===
merged.to_csv(output_file, index=False)
print(f"Filtered and merged results written to {output_file}")
