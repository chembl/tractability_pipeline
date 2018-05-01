import pandas as pd

df = pd.read_csv('b4_mid.csv')

print(len(df.columns))
df = df.groupby('ensembl_gene_id').apply(list, axis=1, as_index=False)

print(df)
