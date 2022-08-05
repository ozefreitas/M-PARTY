from copy import copy
import pandas as pd


sequences_by_cluster_path = "resources/Data/FASTA/CDHIT/"

def read_clustered_seqs(file: str) -> list:
    """Given a Fasta file, processes the sequences returning a list of all the sequences

    Args:
        file (str): Fasta filename

    Returns:
        list: List containing all and only the sequences from the file
    """
    with open(file, "r") as f:
        Lines = [line.strip("\n") for line in f.readlines()]
        Lines.append("$")
        while("" in Lines):
            Lines.remove("")
        lista_seqs = []
        seqs = ""
        for line in Lines[1:]:
            if not line.startswith(">") and not line.startswith("$"):
                seqs += line
            else:
                lista_seqs.append(seqs)
                seqs = ""
    return lista_seqs


def leave_one_out(seq_list: list) -> list:
    """Receiving a list of sequences, yield another list of all sequences but one, iteratively 
    to all sequences

    Args:
        seq_list (list): A list of sequences

    Yields:
        Iterator[list]: List with one less sequence, for all sequences
    """
    for i in range(len(seq_list)):
        lista_seqs_redut = copy(seq_list)
        seq_out = lista_seqs_redut.pop(i)
        yield lista_seqs_redut, seq_out
    

    

clustered_seqs = read_clustered_seqs(sequences_by_cluster_path + "60-65/1.fasta")
for set_prot in leave_one_out(clustered_seqs):
    print(set_prot)