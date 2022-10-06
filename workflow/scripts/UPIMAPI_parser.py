import pandas as pd
import re


def UPIMAPI_parser(filepath):
    UPIMAPI_outfile = pd.read_csv(filepath, sep="\t")
    return UPIMAPI_outfile

def UPIMAPI_iter_per_sim(dataframe):
    """Given a pandas DataFrame, return a dictionary with a list of sequences form the iteration of the sequence similarity between queries and database sequences.

    Args:
        dataframe (DataFrame): A pandas dataframe with diamond documented columns names as header

    Returns:
        dictionary: A dictionary where the keys are intervals of sequence similarity, and values are lists of UniProtKB queries
    """

    # selecionar colunas com perc. identity juntamente com os IDs das sequencias
    # print(dataframe.columns)
    seq_id = dataframe[["qseqid", "sseqid", "pident"]]
    # print(seq_id)

    # retirar os grupos de enzimas com similaridade de 60% a 90% com incrementos de 5%
    target_enzymes = {}
    for perc in range(60, 86, 5):
        chave = str(perc)+"-"+str(perc+5)
        for index, seq in seq_id.iterrows():
            if seq["pident"] >= perc and seq["pident"] < perc+5:
                ident = re.findall("\|.*\|", seq["qseqid"])
                ident = re.sub("\|", "", ident[0])
                if chave not in target_enzymes.keys():
                    target_enzymes[chave] = [ident]
                else:
                    target_enzymes[chave].append(ident)
    return target_enzymes

def save_as_tsv(dic):
    int_df = pd.DataFrame.from_dict(dic, orient="index")
    # int_df.to_csv("C:/Users/jpsfr/OneDrive/Ambiente de Trabalho/TOOL/PDETool/workflow/Data/Tables/UPIMAPI_results_per_sim.tsv", sep="\t")
    int_df.to_csv(snakemake.output[0], sep="\t")


handle = UPIMAPI_parser(snakemake.input[0])
dicionario_identidades = UPIMAPI_iter_per_sim(handle)
save_as_tsv(dicionario_identidades)
