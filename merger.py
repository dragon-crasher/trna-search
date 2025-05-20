import pandas as pd

df1 = pd.read_csv('/mnt/d/bioinformatics/RNAseq_pipeline/diffexp/trf_sequences.csv')
df2 = pd.read_csv('/mnt/d/bioinformatics/RNAseq_pipeline/diffexp/DE_results_healthy_vs_aggressive.csv')

print(df1.head())
cols = df2.columns.tolist()
cols[0] = "tRF sequence"
df2.columns = cols
print(df2.head())

merged_df = pd.merge(df1, df2, on='tRF sequence', how='inner')
fin_df = merged_df.sort_values(by=['linearFoldChange'])
print(fin_df.head())

fin_df.to_csv('/mnt/d/bioinformatics/RNAseq_pipeline/diffexp/CLL_data.csv', index=False)