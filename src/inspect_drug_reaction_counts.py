import pandas as pd

path = "data/processed/faers_case_level/case_reactions.parquet"

df = pd.read_parquet(path)

print(df.head())
print(df.columns)
print(df.shape)
print(df.dtypes)