import pandas as pd
import re
from docker_run import run_command
import os


def concat_hmmsearch_results(folder: str, output: str):
    """Concatenate .TSV files from hmmsearch.

    Args:
        folder (str): Initial folder with the hmmsearch results in .TSV.
        output (str): Output directory for the concatenated files.
    """
    if os.path.exists(output):
        os.remove(output)
    file_list = []
    for i in os.listdir(folder):
        # if os.path.isfile(i):
            file_list.append(folder + "/" + i)
    df_list = []
    for file in file_list:
        df_list.append(pd.read_csv(file, sep = "\t"))
    df = pd.concat(df_list, ignore_index=True)
    df.to_csv(output, sep = "\t")


def BLAST_parser(filepath: str) -> str:
    BLAST_outfile = pd.read_csv(filepath, sep="\t")
    return BLAST_outfile


def BLAST_iter_per_sim(dataframe):
    seq_id = dataframe[["Accession", "Identities(%)"]]
    # seq_id = dataframe[["qseqid", "sseqid", "pident"]]
    # print(seq_id)
    # retirar os grupos de enzimas com similaridade de 60% a 90% com incrementos de 5%
    target_enzymes = {}
    for perc in range(60, 86, 5):
        chave = str(perc)+"-"+str(perc+5)
        for index, seq in seq_id.iterrows():
            if seq["Identities(%)"] >= perc and seq["Identities(%)"] < perc+5:
                if chave not in target_enzymes.keys():
                    target_enzymes[chave] = [seq["Accession"]]
                else:
                    target_enzymes[chave].append(seq["Accession"])
    return target_enzymes


# concat_hmmsearch_results("C:/Users/jpsfr/OneDrive/'Ambiente de Trabalho'/M-PARTY/M-PARTY/resources/Alignments/PET/BLAST/", "C:/Users/jpsfr/OneDrive/Ambiente de Trabalho/M-PARTY/M-PARTY/resources/Alignments/PET/BLAST/test.tsv")
# concat_hmmsearch_results("/home/josepedro/M-PARTY/resources/Alignments/PET/BLAST", "/home/josepedro/M-PARTY/resources/Alignments/PET/BLAST/test.tsv")
# handle = BLAST_parser("/home/josepedro/M-PARTY/resources/Alignments/PET/BLAST/test.tsv")
# dicionario = BLAST_iter_per_sim(handle)
# save_as_tsv(dicionario, "/home/josepedro/M-PARTY/resources/Data/Tables/BLAST_results_per_sim.tsv")
# df = BLAST_parser("/home/josepedro/M-PARTY/resources/Data/Tables/BLAST_results_per_sim.tsv")
# print(df)