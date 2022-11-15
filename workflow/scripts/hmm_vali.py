from copy import copy
import pandas as pd
from docker_run import run_command, docker_run_tcoffee
from hmmsearch_run import run_hmmsearch
from hmm_process import read_hmmsearch_table, get_e_values
import re
import os
import sys
from pathlib import Path
import subprocess
import shutil
sys.path.append("/".join(sys.path[1].split("/")[:-2]))


def make_paths_dic(db_name:str) -> dict:
    """Function that will resolve all folder paths needed for the validation workflow, resolving the
    system path the files are in and including the database name for the models

    Args:
        db_name (str): database name given by the --hmm_db_name option when running PlastEDMA

    Returns:
        dict: a dictionary with an identifier of each path as key and the corresponding path as values
    """
    caminho0 = sys.path[0]
    path_dict = {
        "sequences_by_cluster_path": f'{caminho0}/resources/Data/FASTA/{db_name}/CDHIT/',
        "HMM_directory": f'{caminho0}/resources/Data/HMMs/{db_name}/After_tcoffee_UPI/',
        "vali_directory": f'{caminho0}/resources/Data/FASTA/{db_name}/HMM_oneout/',
        "eliminated_seqs_dir": f'{caminho0}/resources/Data/FASTA/{db_name}/HMM_outseqs/',
        "alignments_test_dir": f'{caminho0}/resources/Alignments/{db_name}/MultipleSequencesAlign/T_Coffee_HMMVal/',
        "hmm_recon_dir": f'{caminho0}/resources/Data/HMMs/{db_name}/HMM_reconstructed/',
        "hmmsearch_results_dir": f'{caminho0}/resources/Data/HMMs/{db_name}/reconstructed_HMMsearch_results/',
        "neg_control_dir": f'{caminho0}/resources/Data/HMMs/{db_name}/negative_cont_HMMsearch_results/',
        "hmmsearch_other_seqs_dir": f'{caminho0}/resources/Data/HMMs/{db_name}/other_seqs_HMMsearch_results/',
        "validated_models_dir": f'{caminho0}/resources/Data/HMMs/{db_name}/validated_HMM/' 
        }
    return path_dict


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


def concat_fasta(hmm_number: str, threshold: str, path_dictionary: dict) -> str:
    """Function to concatenat fasta files, excluding a single file given by the threshold of similarity and 
    the cluster number, that also characterizes the hmm number, for the hmmsearch run of each reconstructed HMM
    against all the other sequences not belonging to that HMM

    Args:
        path_dictionay (dict): dictionary containing the paths to all directories needed for validation
    with the database name given in the --hmm_db_name when running PlastEDMA
        hmm_number (str): number of the built hmm.
        threshold (str): string number of the similarity threshold.
    
    Returns:
        filename (str): string for the created filename. If file already exists, only returns the name of the file,
    not recreating the file all over again
    """
    filename = path_dictionary["hmmsearch_other_seqs_dir"] + threshold + "/" + hmm_number + "_out.fasta"
    fasta_out = path_dictionary["sequences_by_cluster_path"] + threshold + "/" + hmm_number + ".fasta"
    if os.path.exists(filename):
        return filename
    else:
        p = os.listdir(path_dictionary["sequences_by_cluster_path"])
        Path(path_dictionary["hmmsearch_other_seqs_dir"]).mkdir(parents = True, exist_ok = True)
        for thresh in p:
            if not os.path.exists(path_dictionary["hmmsearch_other_seqs_dir"] + thresh):
                os.mkdir(path_dictionary["hmmsearch_other_seqs_dir"] + thresh)
        with open(filename, 'w') as wf:
            for thresh in p:
                path = os.path.join(path_dictionary["sequences_by_cluster_path"], thresh)
                if os.path.isdir(path):
                    for file in file_generator(path, full_path = True):
                        if file == fasta_out:
                            continue
                        with open(file) as f:
                            for line in f:
                                wf.write(line)
                        f.close()
        return filename


def leave_one_out(thresholds: list, path_dictionary: dict):
    """Function that executes all the steps for the HMM validation and filtration with leave-one-out cross 
    validation. Only performs the steps with the removed sequences from the models. Does not return anything,
    only the hmmsearch results.

    Args:
        path_dictionay (dict): dictionary containing the paths to all directories needed for validation
    with the database name given in the --hmm_db_name when running PlastEDMA
        thresholds (list): A list of the thresholds used for similarity separation after UPIMAPI in order to 
    create the correspondent subdirectories
    """
    Path(path_dictionary["vali_directory"]).mkdir(parents = True, exist_ok = True)
    Path(path_dictionary["alignments_test_dir"]).mkdir(parents = True, exist_ok = True)
    Path(path_dictionary["hmm_recon_dir"]).mkdir(parents = True, exist_ok = True)
    Path(path_dictionary["hmmsearch_results_dir"]).mkdir(parents = True, exist_ok = True)
    Path(path_dictionary["eliminated_seqs_dir"]).mkdir(parents = True, exist_ok = True)
    for thresh in thresholds:
        if not os.path.exists(path_dictionary["vali_directory"] + thresh):
            os.mkdir(path_dictionary["vali_directory"] + thresh)
        if not os.path.exists(path_dictionary["alignments_test_dir"] + thresh):
            os.mkdir(path_dictionary["alignments_test_dir"] + thresh)
        if not os.path.exists(path_dictionary["hmm_recon_dir"] + thresh):
            os.mkdir(path_dictionary["hmm_recon_dir"] + thresh)
        if not os.path.exists(path_dictionary["hmmsearch_results_dir"] + thresh):
            os.mkdir(path_dictionary["hmmsearch_results_dir"] + thresh)
        if not os.path.exists(path_dictionary["eliminated_seqs_dir"] + thresh):
            os.mkdir(path_dictionary["eliminated_seqs_dir"] + thresh)
        path = os.path.join(path_dictionary["sequences_by_cluster_path"], thresh)
        if os.path.isdir(path):
            for file in file_generator(path):
                clust_seqs = read_clustered_seqs(path + "/" + file)
                run = 0
                passed = 0
                for set_prot, out_seq in removing_one(clust_seqs):
                    inter =  f'{file.split(".")[0]}_oneless_{run}.fasta'
                    out = f'{file.split(".")[0]}_outseq_{run}.fasta'
                    if os.path.exists(f'{path_dictionary["hmm_recon_dir"] + thresh + "/" + inter.split(".")[0]}.hmm'):
                        run += 1
                        continue
                    write_interfile(path_dictionary["vali_directory"] + thresh + "/" + inter, set_prot)
                    write_interfile(path_dictionary["eliminated_seqs_dir"] + thresh + "/" + out, out_seq, out_sequence = True)
                    try:
                        run_command(f't_coffee`{path_dictionary["vali_directory"] + thresh + "/" + inter}`-output`clustalw_aln`-outfile`{path_dictionary["alignments_test_dir"] + thresh + "/" + inter.split(".")[0]}.clustal_aln`-type`PROTEIN',
                                    sep = "`")
                        # docker_run_tcoffee(f'{sys.path[-1]}/:/data/', 
                        #                     vali_directory + thresh + "/" + inter, 
                        #                     "clustal_aln", 
                        #                     alignments_test_dir + thresh + "/" + f'{inter.split(".")[0]}')
                        run_command(f'hmmbuild {path_dictionary["hmm_recon_dir"] + thresh + "/" + inter.split(".")[0]}.hmm {path_dictionary["alignments_test_dir"] + thresh + "/" + inter.split(".")[0]}.clustal_aln')
                        run_hmmsearch(path_dictionary["eliminated_seqs_dir"] + thresh + "/" + out,
                                        f'{path_dictionary["hmm_recon_dir"] + thresh + "/" + inter.split(".")[0]}.hmm', 
                                        f'{path_dictionary["hmmsearch_results_dir"] + thresh + "/search_" + file.split(".")[0]}_hmm_{run}_seq.tsv', 
                                        out_type = "tsv")
                    except:
                        run += 1
                        print("SOMETHING IS NOT WORKING!!!!!!!!!!")
                        continue
                    df = read_hmmsearch_table(f'{path_dictionary["hmmsearch_results_dir"] + thresh + "/search_" + file.split(".")[0]}_hmm_{run}_seq.tsv')
                    df = check_eval(df)
                    if df is not None:
                        passed += 1
                    run += 1
                # recall divide pelo numero total de sequencias no hmm original (sem leave-one-out)
                hmm_recall = calc_recall(passed, get_number_seqs(f'{path_dictionary["HMM_directory"] + thresh + "/" + file.split(".")[0]}.hmm'))


def negative_control(path_dictionary: dict, database: str = None):
    """Function that executes all the steps for the HMM validation for the negative control
    sequences, in this case, with gut metagenome protein from human. 
    Only performs the steps against the previously built models with one less sequence. Does not return anything, 
    only the hmmsearch results.

    Args:
        path_dictionay (dict): dictionary containing the paths to all directories needed for validation
    with the database name given in the --hmm_db_name when running PlastEDMA
        database (str, optional): A path for a fasta file with a set of protein sequences to serve as negative control.
    Defaults to None. Only for testing purposes.
    """
    if database:
        controlo = database
    else:
        controlo = f'{sys.path[0]}/resources/Data/FASTA/human_gut_metagenome.fasta'
    p = os.listdir(path_dictionary["sequences_by_cluster_path"])
    Path(path_dictionary["neg_control_dir"]).mkdir(parents = True, exist_ok = True)
    for thresh in p:
        if not os.path.exists(path_dictionary["neg_control_dir"] + thresh):
            os.mkdir(path_dictionary["neg_control_dir"] + thresh)
        path = os.path.join(path_dictionary["hmm_recon_dir"], thresh)
        if os.path.isdir(path):
            for hmm in file_generator(path):
                run_hmmsearch(controlo,
                                path + "/" + hmm, 
                                f'{path_dictionary["neg_control_dir"] + thresh}/search_{hmm.split(".")[0]}_{controlo.split(".")[0].split("/")[-1]}.tsv', 
                                out_type = "tsv")


def search_other_seqs(path_dictionary: dict):
    """Function that executes all the steps for the HMM validation with the search against all sequences not included
    in the currently analysing HMM. Does not return anything, only the hmmsearch results. For each threshold full run,
    removes all intermediate fasta files created.
    """
    p = os.listdir(path_dictionary["hmm_recon_dir"])
    for thresh in p:
        path = os.path.join(path_dictionary["hmm_recon_dir"], thresh)
        if os.path.isdir(path):
            for hmm in file_generator(path):
                hmm_num = hmm.split("_")[0]
                target = concat_fasta(hmm_num, thresh, path_dictionary)
                # nome será search_oneless_{sequencia que esta de fora}_clustout{numero do hmm (que é o mesmo do cluster que
                # lhe deu origem e por isso é também a sequencia que ficará de fora)}
                run_hmmsearch(target,
                                path + "/" + hmm, 
                                f'{path_dictionary["hmmsearch_other_seqs_dir"] + thresh}/search_{"_".join(hmm.split("_")[1:]).split(".")[0]}_clustout_{hmm_num}.tsv', 
                                out_type = "tsv")
            for item in os.listdir(path_dictionary["hmmsearch_other_seqs_dir"] + thresh):
                if item.endswith(".fasta"):
                    delete_inter_files(os.path.join(path_dictionary["hmmsearch_other_seqs_dir"] + thresh, item))


def exec_testing(thresholds: list, path_dictionary: dict, database: str = None):
    """Function that executes all the steps for the HMM validation and filtration with leave-one-out cross 
    validation. Does not return anything, just write the final models.

    Args:
        path_dictionay (dict): dictionary containing the paths to all directories needed for validation
    with the database name given in the --hmm_db_name when running PlastEDMA
        thresholds (list): the list of thresholds for the leave-one-out
        database (str, optional): a fasta filename with the sequence database for negative control validation
    """
    leave_one_out(thresholds, path_dictionary)
    negative_control(path_dictionary, database)
    search_other_seqs(path_dictionary)


def hmm_filtration(path_dictionary: dict):
    """Function that will evaluate the results of all three steps of validation - "leave-one-out", "negative
    control" and "search other seqs" and decide which models pass. For this to happen, each reconstructed hmm
    with one less sequence must: all evalues from the recalled sequences must be lower than the minimun evalue 
    from all the other tests (negative control and search agaisnt other sequences).
    Function must differentiate between each hmm in terms of the sequences left out. 

    Args:
        path_dictionay (dict): dictionary containing the paths to all directories needed for validation
        with the database name given in the --hmm_db_name when running PlastEDMA

    Returns:
        false_positives (list): A list containing the number of the models (and respective threshold) which did 
        not passed the validations check.
    """
    # guardar e-values da primeira etapa "leave-one-out"
    p = os.listdir(path_dictionary["hmmsearch_results_dir"])
    # dicionario que irá guardar todos os evalues de cada sequencia recalled para cada hmm com menos uma seq
    eval_per_hmms = {}
    for thresh in p:
        path = os.path.join(path_dictionary["hmmsearch_results_dir"], thresh)
        if os.path.isdir(path):
            for file in file_generator(path):
                # print(path + file)
                hmm_num = file.split("_")[1]
                if f'{thresh}_{hmm_num}' not in eval_per_hmms:
                    eval_per_hmms[f'{thresh}_{hmm_num}'] = []
                df = read_hmmsearch_table(path + "/" + file)
                ev = check_eval(df)
                if ev is not None:
                    eval_per_hmms[f'{thresh}_{hmm_num}'].append(ev)

    # guardar e-values da etapa do controlo negativo
    p = os.listdir(path_dictionary["neg_control_dir"])
    # dicionario que irá guardar todos os evalues de cada sequencia recalled para cada hmm com menos uma seq
    cont_neg_eval_per_hmms = {}
    for thresh in p:
        path = os.path.join(path_dictionary["neg_control_dir"], thresh)
        if os.path.isdir(path):
            for file in file_generator(path):
                # print(path + file)
                hmm_num = file.split("_")[1]
                if f'{thresh}_{hmm_num}' not in cont_neg_eval_per_hmms:
                    cont_neg_eval_per_hmms[f'{thresh}_{hmm_num}'] = []
                df = read_hmmsearch_table(path + "/" + file)
                # print(df)
                ev = check_min_eval(df)
                if ev is not None:
                    cont_neg_eval_per_hmms[f'{thresh}_{hmm_num}'].append(ev)

    # guardar e-values da etapa da procura contra todas as outras seqs
    p = os.listdir(path_dictionary["hmmsearch_other_seqs_dir"])
    # dicionario que irá guardar todos os evalues de cada sequencia recalled para cada hmm com menos uma seq
    other_seqs_eval_per_hmms = {}
    for thresh in p:
        path = os.path.join(path_dictionary["hmmsearch_other_seqs_dir"], thresh)
        if os.path.isdir(path):
            for file in file_generator(path):
                # print(path + file)
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
    # print(eval_per_hmms)
    # print(cont_neg_eval_per_hmms)
    # print(other_seqs_eval_per_hmms)
    for hmm in eval_per_hmms.keys():
        # print(hmm)
        strict_passed = 0
        for eval in range(len(eval_per_hmms[hmm])):
            # print(eval_per_hmms[hmm])
            if float(eval_per_hmms[hmm][eval]) < float(cont_neg_eval_per_hmms[hmm][eval]) and \
            float(eval_per_hmms[hmm][eval]) < float(other_seqs_eval_per_hmms[hmm][eval]):
                strict_passed += 1
            else:
                continue
        # print(hmm)
        hmm_strict_recall = calc_recall(strict_passed, get_number_seqs(f'{path_dictionary["HMM_directory"] + hmm.split("_")[0] + "/" + hmm.split("_")[1]}.hmm'))
        if hmm_strict_recall > 80:
            continue
        else:
            false_positives.append(hmm)
    return false_positives


def remove_fp_models(list_fp: list, path_dictionary: dict):
    """Function that will copy the models checked by validtion to a final folder.

    Args:
        path_dictionay (dict): dictionary containing the paths to all directories needed for validation
        with the database name given in the --hmm_db_name when running PlastEDMA
        list_fp (list): A list containing the number of the models (and respective threshold) which did 
        not passed the validations check.
    """
    Path(path_dictionary["validated_models_dir"]).mkdir(parents = True, exist_ok = True)
    p = os.listdir(path_dictionary["HMM_directory"])
    for thresh in p:
        path = os.path.join(path_dictionary["HMM_directory"], thresh)
        if os.path.isdir(path):
            if not os.path.exists(path_dictionary["validated_models_dir"] + thresh):
                os.mkdir(path_dictionary["validated_models_dir"] + thresh)
            for file in file_generator(path):
                if f'{thresh}_{file.split(".")[0]}' not in list_fp:
                    shutil.copy(path + "/" + file, path_dictionary["validated_models_dir"] + thresh + "/")


def concat_final_model(path_dictionary: dict):
    """Function that will cancatenate the models from each threshold in a single hmm file.

    Args:
        path_dictionay (dict): dictionary containing the paths to all directories needed for validation
        with the database name given in the --hmm_db_name when running PlastEDMA
    """
    p = os.listdir(path_dictionary["validated_models_dir"])
    for thresh in p:
        path = os.path.join(path_dictionary["validated_models_dir"], thresh)
        if os.path.isdir(path):
            with open(path_dictionary["validated_models_dir"] + thresh + ".hmm", "w") as wf:
                for file in file_generator(path):
                    with open(path_dictionary["validated_models_dir"] + thresh + "/" + file, "r") as f:
                        Lines = f.readlines()
                        for line in Lines:
                            wf.write(line)
                    f.close()
            wf.close()


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
# negative_control(database = "resources/Data/FASTA/polymerase_DB.fasta")
# search_other_seqs()
# exec_testing(database = "resources/Data/FASTA/polymerase_DB.fasta")
# a_sair = hmm_filtration()
# print(a_sair)
# remove_fp_models(a_sair)
# concat_final_model()