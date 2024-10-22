import pandas as pd
import anndata as ad

def msdial2anndata(filepath):
    adata = to_anndata(pd.read_csv(filepath, sep='\t', header=None))
    return adata

def to_anndata(df: pd.DataFrame) -> ad.AnnData:
    first_row = df.iloc[0]
    var_start_index = first_row[first_row == 'Class'].index[0]
    df_row = df.iloc[:, 0:var_start_index+1]
    df_col = df.iloc[:, var_start_index+1:]

    df_col_avgstd_removed = df_col.loc[:, ~df_col.iloc[3].str.contains('Average|Stdev')]
    counts = df_col_avgstd_removed.iloc[5:]
    counts_converted = counts.applymap(lambda x: pd.to_numeric(x, errors='coerce'))

    adata = ad.AnnData(counts_converted)
    adata.obs_names = [f"AlignmentID_{i:d}" for i in range(adata.n_obs)]
    adata.var_names = df_col_avgstd_removed.iloc[4]

    # prompt: remove the first 4 rows from df_row
    df_row = df_row.iloc[4:]

    # prompt: set the first row as column names of df_row
    df_row.columns = df_row.iloc[0]
    # prompt: remove the first row from df_row
    df_row = df_row.drop(df_row.index[0])
    # prompt: remove RangeIndex from df_row
    df_row = df_row.reset_index(drop=True)
    # prompt: drop `Alignment ID` column from df_row
    df_row = df_row.drop('Alignment ID', axis=1)

    adata.obs = df_row
    adata.var["Class"] = list(df_col_avgstd_removed.iloc[0])
    adata.var["File type"] = list(df_col_avgstd_removed.iloc[1])
    adata.var["Injection order"] = list(df_col_avgstd_removed.iloc[2])
    adata.var["Batch ID"] = list(df_col_avgstd_removed.iloc[3])

    return adata
