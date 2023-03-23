import pandas as pd
import re
from command_run import run_command


def build_diamond_DB(input_fasta: str, output_path: str, verbose: bool = False) -> str:
	"""Builds a dmnd database file from a fasta file to run with DIAMOND.

	Args:
		input_fasta (str): path for the fasta file.
		output_path (str): path for the output diamond database.
		verbose (bool): prints aditional information. Defaults to False.

	Returns:
		str: name of the resulting file.
	"""
	dmnd_dbname = f'{output_path}/{input_fasta.split("/")[-1].split(".")[0]}'
	if verbose:
		print("Building binary DIAMOND database file...\n")
	run_command(f'diamond`makedb`--in`{input_fasta}`-d`{dmnd_dbname}', sep = "`")
	if verbose:
		print("Done\n")
	return dmnd_dbname + ".dmnd"


def run_DIAMOND(query: str, outpath: str, database: str, threads: int) -> str:
    """Function to run DIAMOND by bash with the given arguments from M-PARTY.

    Args:
        query (str): Database path to be searched on by DIAMOND.
        outpath (str): Output directory.
        database (str): Sequence list to be searched against the database.
        threads (int): Number of threads.

    Returns:
        str: Path to the final TSV file.
    """
    run_command(f'diamond`blastp`-q`{query}`-o`{outpath}`-d`{database}`--threads`{threads}`--very-sensitive`--outfmt`6`--unal`1`-b`0.36036930084228513`-c`4`--evalue`0.001', sep = "`")
    return outpath


def DIAMOND_parser(filepath: str):
    DIAMOND_outfile = pd.read_csv(filepath, sep="\t")
    DIAMOND_outfile.columns = ["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore"]
    return DIAMOND_outfile


def DIAMOND_iter_per_sim(dataframe: pd.DataFrame) -> dict:
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

