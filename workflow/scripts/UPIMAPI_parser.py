import pandas as pd
import re
from command_run import run_command


def run_UPIMAPI(query: str, outpath: str, upi_database: str, threads: int) -> str:
    """Function to run UPIMAPI by bash with the given arguments from M-PARTY.

    Args:
        query (str): Database path to be searched on by DIAMOND.
        outpath (str): Output directory name with filename.
        upi_database (str): Sequence list to be searched against the database.
        threads (int): Number of threads.

    Returns:
        str: Path to the final TSV file.
    """
    run_command(f'upimapi.py`-i`{query}`-o`{outpath}`--database`{upi_database}`-t`{threads}', sep = "`")
    return outpath + "/UPIMAPI_results.tsv"


def UPIMAPI_parser(filepath: str):
    UPIMAPI_outfile = pd.read_csv(filepath, sep="\t")
    return UPIMAPI_outfile


def UPIMAPI_iter_per_sim(dataframe: pd.DataFrame) -> dict:
    """Given a pandas DataFrame, return a dictionary with a list of sequences form the iteration of the sequence similarity between queries and database sequences.

    Args:
        dataframe (pd.DataFrame): A pandas dataframe with diamond documented columns names as header.

    Returns:
        dict: A dictionary where the keys are intervals of sequence similarity, and values are lists of UniProtKB queries.
    """
    # selecionar colunas com perc. identity juntamente com os IDs das sequencias
    # print(dataframe.columns)
    seq_id = dataframe[["qseqid", "sseqid", "pident"]]
    # print(seq_id)
    # retirar os grupos de enzimas com similaridade de no minimo 60% a 90% com incrementos de 5%
    target_enzymes = {}
    for perc in range(60, 91, 5):
        # chave = str(perc)+"-"+str(perc+5)
        for index, seq in seq_id.iterrows():
            # print(type(seq["pident"]))
            # if seq["pident"] >= perc and seq["pident"] < perc+5:
            if seq["pident"] >= perc:
                ident = re.findall("\|.*\|", seq["qseqid"])
                ident = re.sub("\|", "", ident[0])
                if perc not in target_enzymes.keys():
                    target_enzymes[perc] = [ident]
                else:
                    target_enzymes[perc].append(ident)
    return target_enzymes


# handle = UPIMAPI_parser(snakemake.input[0])
# dicionario_identidades = UPIMAPI_iter_per_sim(handle)
# save_as_tsv(dicionario_identidades)
