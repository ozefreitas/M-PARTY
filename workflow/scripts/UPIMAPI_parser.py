import pandas as pd
import re
from docker_run import run_command


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
    # run_command(f'upimapi.py`-i`{query}`-o`{outdir}`--database`{upi_database}`-t`{threads}', sep = "`")
    run_command(f'diamond`blastp`-q`{query}`-o`{outpath}`-d`{upi_database}`--threads`{threads}`--very-sensitive`--outfmt`6`--unal`1`--max-target-seqs`1`-b`0.36036930084228513`-c`4`--evalue`0.001', sep = "`")
    # return outdir + "/UPIMAPI_results.tsv"
    return outpath

def UPIMAPI_parser(filepath: str):
    UPIMAPI_outfile = pd.read_csv(filepath, sep="\t")
    UPIMAPI_outfile.columns = ["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore"]
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
    print(seq_id)
    # retirar os grupos de enzimas com similaridade de 60% a 90% com incrementos de 5%
    target_enzymes = {}
    for perc in range(60, 86, 5):
        chave = str(perc)+"-"+str(perc+5)
        for index, seq in seq_id.iterrows():
            print(type(seq["pident"]))
            if seq["pident"] >= perc and seq["pident"] < perc+5:
                ident = re.findall("\|.*\|", seq["qseqid"])
                ident = re.sub("\|", "", ident[0])
                if chave not in target_enzymes.keys():
                    target_enzymes[chave] = [ident]
                else:
                    target_enzymes[chave].append(ident)
    return target_enzymes

def save_as_tsv(dic: dict, out_path: str):
    int_df = pd.DataFrame.from_dict(dic, orient="index")
    int_df.to_csv(out_path, sep="\t")


# handle = UPIMAPI_parser(snakemake.input[0])
# dicionario_identidades = UPIMAPI_iter_per_sim(handle)
# save_as_tsv(dicionario_identidades)
