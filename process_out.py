import pandas as pd

df = pd.read_csv('out_buckets2.csv')
print(df)


df.to_csv('processed_buckets.csv')