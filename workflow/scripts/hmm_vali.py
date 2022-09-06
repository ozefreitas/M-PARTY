from copy import copy
import pandas as pd
from docker_run import run_command, docker_run_tcoffee
from hmmsearch_run import run_hmmsearch
from hmm_process import read_hmmsearch_table, get_e_values
import re
import os
import sys
import subprocess
sys.path.append("/".join(sys.path[0].split("/")[:-2]))


# passar para dicionario!!
sequences_by_cluster_path = "resources/Data/FASTA/CDHIT/"
HMM_directory = "resources/Data/HMMs/After_tcoffee_UPI/"
vali_directory = "resources/Data/FASTA/HMM_oneout/"
eliminated_seqs_dir = "resources/Data/FASTA/HMM_outseqs/"
alignments_test_dir = "resources/Alignments/MultipleSequencesAlign/T_Coffee_HMMVal/"
hmm_recon_dir = "resources/Data/HMMs/HMM_reconstructed/"
hmmsearch_results_dir = "resources/Data/HMMs/reconstructed_HMMsearch_results/"
neg_control_dir = "resources/Data/HMMs/negative_cont_HMMsearch_results/"
hmmsearch_other_seqs_dir = "resources/Data/HMMs/other_seqs_HMMsearch_results/"


def file_generator(path: str, full_path: bool = False) -> str:
    """Function that yield the name of all and only files inside a directory in the given path, for iteration purposes.

    Args:
        path (str): Path for the folder to be analyzed.
        full_path (bool, optional): If wants to return the full path of the file. Defaults to False.

    Yields:
        str: Yield the name of each file inside the given directory.
    """
    for file in os.listdir(path):
        if os.path.isfile(os.path.join(path, file)):
            if full_path:
                yield os.path.join(path, file)
            else:
                yield file


def delete_inter_files(file_path: str):
    """Quick function to delete intermediate files.

    Args:
        file_path (str): relative path to the file.

    Raises:
        FileNotFoundError: raised if file does not exist.
    """
    if os.path.exists(file_path):
        os.remove(file_path)
    else:
        raise FileNotFoundError(file_path)


def read_clustered_seqs(file: str) -> list:
    """Given a Fasta file, processes the sequences returning a list of all the sequences.

    Args:
        file (str): Fasta filename.

    Returns:
        list: List containing all the sequences from the file, a full sequence (with ID) by element.
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


def removing_one(seq_list: list) -> list:
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
    

def write_interfile(filename: str, seq_list: list, out_sequence: bool = False):
    """Will write the intermediate fasta files with one less sequence, so they can then be aligned,
    and build the HMM from them.

    Args:
        filename (str): Name of the output intermediate file (will be autommatically generated).
        seq_list (list): List of sequences from leave-one-out function.
        out_sequence (bool, optional): If the file to be writen is the sequences that was leaved out. Defaults to False.
    """
    if out_sequence:
        with open(filename, "w") as wf:
            search = seq_list.index("\n")
            wf.write(seq_list[:search])
            fasta_form = []
            for i in range(0, len(seq_list) - search, 60):
                fasta_form.append(seq_list[i + search: i + search + 60] + "\n")
            fasta_form.append("\n")
            for x in fasta_form:
                wf.write(x)
        wf.close()
    else:
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


def get_number_seqs(model: str, file_type: str = "HMM") -> int:
    """Given a file with HMM related info, can be a fasta that will construct the HMM of the HMM itself, returns
    the number of sequences composing that model.

    Args:
        model (str): model filename.
        file_type (str, optional): type of file between HMM and FASTA. Defaults to "HMM".

    Returns:
        int: number of sequences in the model.
    """
    with open(model, "r") as wf:
        Lines = []
        for _ in range(16):
            Lines.append(wf.readline())
    wf.close()
    for line in Lines:
        number = re.findall("NSEQ.*", line)
        if number != []:
            number = re.findall("[0-9]", number[0])
            break
    return int("".join(number))


def check_eval(dataframe: pd.DataFrame) -> int:
    """Given a Dataframe from the conversion of the hmmsearch results, reads the single e-value from that
    dataframe and decides wether is acceptable to proced. Returns the e-value if e-value threshold 
    is higher than the defined and returns None if otherwise.

    Args:
        dataframe (pd.DataFrame): A Dataframe from the conversion of the hmmsearch results txt file to pandas.
    
    Returns:
        int: The e-value
    """
    eval = get_e_values(dataframe, to_list = True)[0]
    if int(float(eval)) < 10:
        return eval
    else:
        return None


def check_min_eval(dataframe: pd.DataFrame) -> int:
    """Given a Dataframe from the conversion of the hmmsearch results, reads the e-values from that
    dataframe and decides wether is acceptable to proced. Returns the minimum e-value from the list.
    If e-value from the tsv is not available because no value was found, returns an abnormally large number.

    Args:
        dataframe (pd.DataFrame): A Dataframe from the conversion of the hmmsearch results txt file to pandas.

    Returns:
        int: The minimum e-value.
    """
    eval = get_e_values(dataframe, to_list = True)
    try:
        return min(eval)
    except:
        return 1


def calc_recall(recalled: int, hmm_size: int) -> int:
    return (recalled/hmm_size) * 100


def calc_strict_recall(strictly_recalled: int, hmm_size: int) -> int:
    return (strictly_recalled/hmm_size) * 100


def concat_fasta(hmm_number: str, threshold: str) -> str:
    """Function to concatenat fasta files, excluding a single file given by the threshold of similarity and 
    the cluster number, that also characterizes the hmm number, for the hmmsearch run of each reconstructed HMM
    against all the other sequences not belonging to that HMM

    Args:
        hmm_number (str): number of the built hmm.
        threshold (str): string number of the similarity threshold.
    
    Returns:
        filename (str): string for the created filename. If file already exists, only returns the name of the file,
        not recreating the file all over again
    """
    filename = hmmsearch_other_seqs_dir + threshold + "/" + hmm_number + "_out.fasta"
    fasta_out = sequences_by_cluster_path + threshold + "/" + hmm_number + ".fasta"
    if os.path.exists(filename):
        return filename
    else:
        p = os.listdir(sequences_by_cluster_path)
        for thresh in p:
            if not os.path.exists(hmmsearch_other_seqs_dir):
                os.mkdir(hmmsearch_other_seqs_dir)
            if not os.path.exists(hmmsearch_other_seqs_dir + thresh):
                os.mkdir(hmmsearch_other_seqs_dir + thresh)
        with open(filename, 'w') as wf:
            for thresh in p:
                path = os.path.join(sequences_by_cluster_path, thresh)
                if os.path.isdir(path):
                    for file in file_generator(path, full_path = True):
                        if file == fasta_out:
                            continue
                        with open(file) as f:
                            for line in f:
                                wf.write(line)
                        f.close()
        return filename


def leave_one_out():
    """Function that executes all the steps for the HMM validation and filtration with leave-one-out cross 
    validation. Only performs the steps with the removed sequences from the models. Does not return anything,
    only the hmmsearch results.
    """
    p = os.listdir(sequences_by_cluster_path)
    for thresh in p:
        if not os.path.exists(vali_directory + thresh):
            os.mkdir(vali_directory + thresh)
        if not os.path.exists(alignments_test_dir + thresh):
            os.mkdir(alignments_test_dir + thresh)
        if not os.path.exists(hmm_recon_dir + thresh):
            os.mkdir(hmm_recon_dir + thresh)
        if not os.path.exists(hmmsearch_results_dir + thresh):
            os.mkdir(hmmsearch_results_dir + thresh)
        if not os.path.exists(eliminated_seqs_dir + thresh):
            os.mkdir(eliminated_seqs_dir + thresh)
        path = os.path.join(sequences_by_cluster_path, thresh)
        if os.path.isdir(path):
            for file in file_generator(path):
                clust_seqs = read_clustered_seqs(path + "/" + file)
                run = 0
                passed = 0
                for set_prot, out_seq in removing_one(clust_seqs):
                    inter =  f'{file.split(".")[0]}_oneless_{run}.fasta'
                    out = f'{file.split(".")[0]}_outseq_{run}.fasta'
                    if os.path.exists(f'{hmm_recon_dir + thresh + "/" + inter.split(".")[0]}.hmm'):
                        run += 1
                        continue
                    write_interfile(vali_directory + thresh + "/" + inter, set_prot)
                    write_interfile(eliminated_seqs_dir + thresh + "/" + out, out_seq, out_sequence = True)
                    try:
                        docker_run_tcoffee(f'{sys.path[-1]}/:/data/', 
                                            vali_directory + thresh + "/" + inter, 
                                            "clustal_aln", 
                                            alignments_test_dir + thresh + "/" + f'{inter.split(".")[0]}')
                        run_command(f'hmmbuild {hmm_recon_dir + thresh + "/" + inter.split(".")[0]}.hmm {alignments_test_dir + thresh + "/" + inter.split(".")[0]}.clustal_aln')
                        run_hmmsearch(eliminated_seqs_dir + thresh + "/" + out,
                                        f'{hmm_recon_dir + thresh + "/" + inter.split(".")[0]}.hmm', 
                                        f'{hmmsearch_results_dir + thresh + "/search_" + file.split(".")[0]}_hmm_{run}_seq.tsv', 
                                        out_type = "tsv")
                    except:
                        run += 1
                        continue
                    df = read_hmmsearch_table(f'{hmmsearch_results_dir + thresh + "/search_" + file.split(".")[0]}_hmm_{run}_seq.tsv')
                    df = check_eval(df)
                    if df is not None:
                        passed += 1
                    run += 1
                # recall divide pelo numero total de sequencias no hmm original (sem leave-one-out)
                hmm_recall = calc_recall(passed, get_number_seqs(f'{HMM_directory + thresh + "/" + file.split(".")[0]}.hmm'))


def negative_control():
    """Function that executes all the steps for the HMM validation for the negative control
    sequences, in this case, with gut metagenome protein from human. 
    Only performs the steps against the previously built models with one less sequence. Does not return anything, 
    only the hmmsearch results.
    """
    controlo = "resources/Data/FASTA/human_gut_metagenome.fasta"
    p = os.listdir(sequences_by_cluster_path)
    for thresh in p:
        if not os.path.exists(neg_control_dir + thresh):
            os.mkdir(neg_control_dir + thresh)
        path = os.path.join(hmm_recon_dir + thresh)
        if os.path.isdir(path):
            for hmm in file_generator(path):
                run_hmmsearch(controlo,
                                path + "/" + hmm, 
                                f'{neg_control_dir + thresh}/search_{hmm.split(".")[0]}_human_gut_metagenome.tsv', 
                                out_type = "tsv")


def search_other_seqs():
    """Function that executes all the steps for the HMM validation with the search against all sequences not included
    in the currently analysing HMM. Does not return anything, only the hmmsearch results. For each threshold full run,
    removes all intermediate fasta files created.
    """
    p = os.listdir(hmm_recon_dir)
    for thresh in p:
        path = os.path.join(hmm_recon_dir + thresh)
        if os.path.isdir(path):
            for hmm in file_generator(path):
                hmm_num = hmm.split("_")[0]
                target = concat_fasta(hmm_num, thresh)
                # nome será search_oneless_{sequencia que esta de fora}_clustout{numero do hmm (que é o mesmo do cluster que
                # lhe deu origem e por isso é também a sequencia que ficará de fora)}
                run_hmmsearch(target,
                                path + "/" + hmm, 
                                f'{hmmsearch_other_seqs_dir + thresh}/search_{"_".join(hmm.split("_")[1:]).split(".")[0]}_clustout_{hmm_num}.tsv', 
                                out_type = "tsv")
            for item in os.listdir(hmmsearch_other_seqs_dir + thresh):
                if item.endswith(".fasta"):
                    delete_inter_files(os.path.join(hmmsearch_other_seqs_dir + thresh, item))


def exec_testing():
    """Function that executes all the steps for the HMM validation and filtration with leave-one-out cross 
    validation. Does not return anything, just write the final models.
    """
    leave_one_out()
    negative_control()
    search_other_seqs()


def hmm_filtration():
    """Function that will evaluate the results of all three steps of validation - "leave-one-out", "negative
    control" and "search other seqs" and decide which models pass. For this to happen, each reconstructed hmm
    with one less sequence must: all evalues from the recalled sequences must be lower than the minimun evalue 
    from all the other tests (negative control and search agaisnt other sequences).
    Function must differentiate between each hmm in terms of the sequences left out. 

    Returns:
        false_positives (list): A list containing the number of the models (and respective threshold) which did 
        not passed the validations check.
    """
    # guardar e-values da primeira etapa "leave-one-out"
    p = os.listdir(hmmsearch_results_dir)
    # dicionario que irá guardar todos os evalues de cada sequencia recalled para cada hmm com menos uma seq
    eval_per_hmms = {}
    for thresh in p:
        path = os.path.join(hmmsearch_results_dir + thresh)
        if os.path.isdir(path):
            for file in file_generator(path):
                hmm_num = file.split("_")[1]
                if f'{thresh}_{hmm_num}' not in eval_per_hmms:
                    eval_per_hmms[f'{thresh}_{hmm_num}'] = []
                df = read_hmmsearch_table(path + "/" + file)
                ev = check_eval(df)
                if ev is not None:
                    eval_per_hmms[f'{thresh}_{hmm_num}'].append(ev)

    # guardar e-values da etapa do controlo neagtivo
    p = os.listdir(neg_control_dir)
    # dicionario que irá guardar todos os evalues de cada sequencia recalled para cada hmm com menos uma seq
    cont_neg_eval_per_hmms = {}
    for thresh in p:
        path = os.path.join(neg_control_dir + thresh)
        if os.path.isdir(path):
            for file in file_generator(path):
                hmm_num = file.split("_")[1]
                if f'{thresh}_{hmm_num}' not in cont_neg_eval_per_hmms:
                    cont_neg_eval_per_hmms[f'{thresh}_{hmm_num}'] = []
                df = read_hmmsearch_table(path + "/" + file)
                ev = check_min_eval(df)
                if ev is not None:
                    cont_neg_eval_per_hmms[f'{thresh}_{hmm_num}'].append(ev)

    # guardar e-values da etape da procura contra todas as outras seqs
    p = os.listdir(hmmsearch_other_seqs_dir)
    # dicionario que irá guardar todos os evalues de cada sequencia recalled para cada hmm com menos uma seq
    other_seqs_eval_per_hmms = {}
    for thresh in p:
        path = os.path.join(hmmsearch_other_seqs_dir + thresh)
        if os.path.isdir(path):
            for file in file_generator(path):
                hmm_num = file.split("_")[-1].split(".")[0]
                if f'{thresh}_{hmm_num}' not in other_seqs_eval_per_hmms:
                    other_seqs_eval_per_hmms[f'{thresh}_{hmm_num}'] = []
                df = read_hmmsearch_table(path + "/" + file)
                ev = check_min_eval(df)
                if ev is not None:
                    other_seqs_eval_per_hmms[f'{thresh}_{hmm_num}'].append(ev)

    # lista de modelos que não passaram a validação
    false_positives = []
    # comparar e-values da primeira etapa com os correspondestes das etapas de validação
    for hmm in eval_per_hmms.keys():
        strict_passed = 0
        for eval in eval_per_hmms[hmm]:
            if eval < cont_neg_eval_per_hmms[hmm][0] and eval < other_seqs_eval_per_hmms[hmm][0]:
                strict_passed += 1
            else:
                continue
        hmm_strict_recall = calc_recall(strict_passed, get_number_seqs(f'{HMM_directory + thresh + "/" + hmm.split("_")[1]}.hmm'))
        if hmm_strict_recall > 80:
            continue
        else:
            false_positives.append(hmm)
    return false_positives


# clustered_seqs = read_clustered_seqs(sequences_by_cluster_path + "60-65/1.fasta")
# # # print(clustered_seqs)
# run = 0
# for set_prot, out_seq in leave_one_out(clustered_seqs):
# #     print(set_prot, out_seq)
#     file = f'qualquercoisa_oneless_{run}'
#     write_interfile(f'{file}.fasta', set_prot)
#     run += 1
# # # print(get_number_seqs(HMM_directory + "60-65/1.hmm"))
#     docker_run_tcoffee(f'{sys.path[-1]}/:/data/', f'{file}.fasta', "clustal_aln", file)
#     run_command(f'hmmbuild {file}.hmm {file}.clustal_aln')
#     run_hmmsearch(out_seq, f'{file}.hmm', f'{file}.tsv', out_type="tsv")
# ola = concat_fasta("1", "60-65")
# delete_inter_files("resources/Data/HMMs/other_seqs_HMMsearch_results/60-65/1out.fasta")

# leave_one_out()
# negative_control()
# search_other_seqs()
# exec()
# hmm_filtration()