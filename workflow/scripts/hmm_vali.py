from copy import copy
import pandas as pd


sequences_by_cluster_path = "resources/Data/FASTA/CDHIT/"

def read_clustered_seqs(file: str) -> list:
    """_summary_

    Args:
        file (str): _description_

    Returns:
        list: _description_
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
    """_summary_

    Args:
        seq_list (list): _description_

    Returns:
        list: _description_

    Yields:
        Iterator[list]: _description_
    """
    for i in range(len(seq_list)):
        lista_seqs_redut = copy(seq_list)
        seq_out = lista_seqs_redut.pop(i)
        yield lista_seqs_redut, seq_out
    

clustered_seqs = read_clustered_seqs(sequences_by_cluster_path + "60-65/1.fasta")
for set_prot in leave_one_out(clustered_seqs):
    print(set_prot)