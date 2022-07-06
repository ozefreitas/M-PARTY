import re
import pandas as pd


def cdhit_parser(txtfile):
    """Receives a text file with a similar format as a FASTA file, and returns a dictionary with the number of the cluster as key and the UniProt ID's for the sequences inside each cluster as value.

    Args:
        txtfile (string): The name for the cluster file generated by CD-HIT

    Returns:
        dictionary: A dictionary with the number of the cluster as key and the UniProt ID's for the sequences inside each cluster as value.
    """
    file = open(txtfile, "r")
    cluster = 0
    seqs_by_cluster = {}
    Lines = file.readlines()
    for line in Lines:
        if line[0] == ">":
            cluster += 1
            seqs_by_cluster[cluster] = []
        else:
            target_seq = re.findall("\|.*\|", line)
            # print(target_seq)
            clean = re.sub("\|", "", target_seq[0])
            seqs_by_cluster[cluster].append(clean)
    return seqs_by_cluster

def counter(clstr_lst, remove_single=True, remove_duplicates=False, tsv_ready=False):
    """Functions receives a dictionary with keys as the number of the cluster and UniProt sequences IDs as values and returns another dictionary
    with the number of sequences per cluster, with the option of removing single sequence clusters.

    Args:
        clstr_lst (dictionary): A dictionary with the number of the cluster as key and the UniProt ID's for the sequences inside each cluster as value.
        remove_single (bool, optional): Decides to remove single sequence clusters. Defaults to True.
        tsv_ready (bool, optional): Only works if remove_single is set to True, and makes a dicitionary of lists, ready to be saved as tsv. Defaults to False.
        remove_duplicates (bool, optional): Decides to remove duplicates from the clusters. Defaults to False.
    Returns:
        dictionary: A dictionary with the number of the cluster as key and the UniProt ID's for the sequences inside each cluster, as well as the size
        of that same cluster as value, in the form of tuple.
    """
    number_seqs_by_cluster = {}
    for k, v in clstr_lst.items():
        if remove_single:
            if len(v) > 1:
                if tsv_ready:
                    number_seqs_by_cluster[k] = v
                else:
                    number_seqs_by_cluster[k] = (v, len(v))
        else:
            number_seqs_by_cluster[k] = (v, len(v))
    if remove_duplicates:
        for k, v in number_seqs_by_cluster.items():
            number_seqs_by_cluster[k] = list(set(v))
    return number_seqs_by_cluster

def save_as_tsv(dic, out_path):
    int_df = pd.DataFrame.from_dict(dic, orient="index")
    int_df.to_csv(out_path, sep="\t")


handle = cdhit_parser(snakemake.input[0])
handle2 = counter(handle, tsv_ready=True, remove_duplicates=True)
save_as_tsv(handle2, snakemake.output[0])
