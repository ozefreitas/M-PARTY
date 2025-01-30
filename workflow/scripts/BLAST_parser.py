import pandas as pd
import re
from workflow.scripts.command_run import run_command
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


def build_blast_DB(input_fasta: str, output_path: str, input_type: str, verbose: bool = False) -> str:
    """Builds a blast database file from a fasta file to run with DIAMOND.

    Args:
        input_fasta (str): path for the fasta file.
        output_path (str): path for the output diamond database.
        input_type (str): Sequence type in file: protein ou nucleic.
        verbose (bool): prints aditional information. Defaults to False.

    Returns:
        str: name of the resulting file.
    """
    blast_dbname = f'{output_path}/{input_fasta.split("/")[-1].split(".")[0]}'
    if input_type in  ["protein", "proteic", "aa", "p", "prot"]:
        input_type = "prot"
    elif input_type in ["nucleic", "nucleotid", "n", "nucl", "nuc"]:
        input_type = "nucl"
    if verbose:
        print("Building BLAST database file...\n")
    run_command(f'makeblastdb`-in`{input_fasta}`-out`{blast_dbname}`-dbtype`{input_type}`-title`BLAST_run`-parse_seqids', sep = "`")
    if verbose:
        print("Done\n")
    return blast_dbname


def run_BLAST(query: str, outpath: str, database: str, threads: int) -> str:
    """Function to run BLAST by bash with the given arguments from M-PARTY.

    Args:
        query (str): Sequence list to be searched against the database.
        outpath (str): Output directory name with filename.
        database (str): Database path to be searched against query sequences.
        threads (str): Number of threads.

    Returns:
        str: Path to the final .TSV file. 
    """
    run_command(f'blastp`-query`{query}`-out`{outpath}`-db`{database}`-num_threads`{threads}`-outfmt`6', sep = "`")
    return outpath


def BLAST_parser(filepath: str) -> str:
    BLAST_outfile = pd.read_csv(filepath, sep="\t")
    BLAST_outfile.columns = ["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore"]
    return BLAST_outfile


def BLAST_iter_per_sim(dataframe):
    # seq_id = dataframe[["Accession", "Identities(%)"]]
    seq_id = dataframe[["qseqid", "sseqid", "pident"]]
    # retirar os grupos de enzimas com similaridade minima de 60% a 90% com incrementos de 5%
    target_enzymes = {}
    for perc in range(60, 91, 5):
        # chave = str(perc)+"-"+str(perc+5)
        for index, seq in seq_id.iterrows():
            # if seq["Identities(%)"] >= perc and seq["Identities(%)"] < perc+5:
            if seq["pident"] >= perc:
                if perc not in target_enzymes.keys():
                    target_enzymes[perc] = [seq["sseqid"]]
                else:
                    target_enzymes[perc].append(seq["sseqid"])
    return target_enzymes


# concat_hmmsearch_results("C:/Users/jpsfr/OneDrive/'Ambiente de Trabalho'/M-PARTY/M-PARTY/resources/Alignments/PET/BLAST/", "C:/Users/jpsfr/OneDrive/Ambiente de Trabalho/M-PARTY/M-PARTY/resources/Alignments/PET/BLAST/test.tsv")
# concat_hmmsearch_results("/home/josepedro/M-PARTY/resources/Alignments/PET/BLAST", "/home/josepedro/M-PARTY/resources/Alignments/PET/BLAST/test.tsv")
# handle = BLAST_parser("/home/josepedro/M-PARTY/resources/Alignments/PET/BLAST/test.tsv")
# dicionario = BLAST_iter_per_sim(handle)
# save_as_tsv(dicionario, "/home/josepedro/M-PARTY/resources/Data/Tables/BLAST_results_per_sim.tsv")
# df = BLAST_parser("/home/josepedro/M-PARTY/resources/Data/Tables/BLAST_results_per_sim.tsv")
# print(df)