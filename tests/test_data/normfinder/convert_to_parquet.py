import polars as pl

df = pl.read_csv("/home/olivier/repositories/nf-core-stableexpression/tests/all_counts.normfinder.csv")
df.write_parquet("/home/olivier/repositories/nf-core-stableexpression/tests/all_counts.normalised.parquet")
