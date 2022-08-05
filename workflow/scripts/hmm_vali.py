from copy import copy
from genericpath import isdir
import pandas as pd
from docker_run import run_command, docker_run_tcoffee
import re
import os
import sys


sequences_by_cluster_path = "resources/Data/FASTA/CDHIT/"
vali_directory = "resources/Data/HMMs/HMM_validation_test/"


def file_generator(path: str, full_path: bool = False) -> str:
    """Function that yield the name of all and only files inside a directory in the given path, for iteration purposes
    Args:
        path (str): Path for the folder to be analyzed

    Yields:
        str: Yield the name of each file inside the given directory
    """

    for file in os.listdir(path):
        if os.path.isfile(os.path.join(path, file)):
            if full_path:
                yield os.path.join(path, file)
            else:
                yield file


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
        seqs = Lines[0]+"\n"
        for line in Lines[1:]:
            if not line.startswith(">") and not line.startswith("$"):
                seqs += line
            else:
                lista_seqs.append(seqs)
                seqs = line+"\n"
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
    

def write_interfile(filename: str, seq_list: list):
    """Will write the intermediate fasta files with one less sequence, so they can then be aligned,
    and build the HMM from them.

    Args:
        filename (str): Name of the output intermediate file (will be autommatically generated)
        seq_list (list): List of sequences from leave-one-out function
    """
    with open(filename, "w") as wf:
        for seq in seq_list:
            search = seq.index("\n")
            wf.write(seq[:search])
            fasta_form = []
            for i in range(0, len(seq) - search, 60):
                fasta_form.append(seq[i + search: i + search + 60] + "\n")
            fasta_form.append("\n")
            for x in fasta_form:
                wf.write(x)
    wf.close()


def exec():
    p = os.listdir(sequences_by_cluster_path)
    for i in p:
        if os.path.isdir(i):
            for file in file_generator(i):
                clust_seqs = read_clustered_seqs(file)
                run = 0
                for set_prot, out_seq in leave_one_out(clust_seqs):
                    file =  f'{file.split(".")[0]}_oneless_{run}.fasta'
                    write_interfile(file, set_prot)
                    docker_run_tcoffee(sys.path[0], file, "clustal_aln", f'{file.split(".")[:-1]}.clustal_aln')
                    run_command(f'hmmbuild {file.split(".")[:-1]}.hmm {file}')
                    run_command(f'hmmsearch ')
                    run += 1


# clustered_seqs = read_clustered_seqs(sequences_by_cluster_path + "60-65/1.fasta")
# # print(clustered_seqs)
# for set_prot, out_seq in leave_one_out(clustered_seqs):
#     print(set_prot, out_seq)
# # write_interfile("qualquercoisa.fasta", clustered_seqs)