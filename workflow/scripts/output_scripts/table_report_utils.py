import pandas as pd

def write_excel(outdir: str, dataframe: pd.DataFrame, ids_per_model: dict) -> None:
    writer = pd.ExcelWriter(outdir + "report_table.xlsx", engine = "openpyxl")
    dataframe.to_excel(writer, sheet_name = "Table_Report", index = 0)
    df1 = pd.DataFrame.from_dict(ids_per_model, orient = "index")
    df1.to_excel(writer, sheet_name = "Model_Sequences")
    writer.save()  # FutureWarning: save is not part of the public API, usage can give unexpected results and will be removed in a future version
    writer.close()

def write_tsv(outdir, table_name, dataframe) -> None:
    dataframe.to_csv(outdir + table_name, sep = "\t")

def write_csv(outdir, table_name, dataframe) -> None:
    dataframe.to_csv(outdir + table_name)

def check_output(type: str, outdir: str, table_name: str, dataframe: pd.DataFrame, ids_per_model: dict) -> None:
    if type == "tsv":
        write_tsv(outdir, table_name, dataframe)
    elif type == "csv":
        write_csv(outdir, table_name, dataframe)
    elif type == "excel":
        write_excel(outdir, dataframe, ids_per_model)
    else:
        raise TypeError(f'Specified table format {type} is not available. Read documentation for --output_type')
