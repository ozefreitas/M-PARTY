#!/usr/bin/env python
# run tool main script without indicating python

"""
M-PARTY - Mining Protein dAtasets for Target EnzYmes

by José Freitas

Oct 2022
"""

import argparse
import sys
sys.path.insert(0, f'{"/".join(sys.path[0].split("/")[:-1])}/share')
sys.path.append(f'{sys.path[1]}/workflow/scripts')
# sys.path.append(f'{sys.path[0]}/M-PARTY')
# print(sys.path)
import os
from pathlib import Path, PureWindowsPath
import time
import yaml
import json
import re
import pandas as pd
from collections import Counter
import glob
import snakemake

from hmmsearch_run import run_hmmsearch
from hmm_process import *
from hmm_vali import concat_final_model, file_generator, exec_testing, hmm_filtration, remove_fp_models, make_paths_dic, delete_inter_files
from UPIMAPI_parser import *
from CDHIT_parser import *
from mparty_util import build_UPI_query_DB, threshold2clusters, get_tsv_files, build_diamond_DB
from BLAST_parser import *


version = "0.2.3"

parser = argparse.ArgumentParser(description="M-PARTY's main script")
parser.add_argument("-i", "--input", help = "input FASTA file containing\
                    a list of protein sequences to be analysed")
parser.add_argument("--input_seqs_db_const", help = "input a FASTA file with a set of sequences from which the user \
                    wants to create the HMM database from scratch")
parser.add_argument("-db", "--database", help = "FASTA database to run against the also user inputted sequences. DIAMOND \
                    is performed in order to expand the data and build the models. M-PARTY has no in-built database for this \
                    matter. If flag is given, download of the default database will start and model built from that. Defaults to UniProt DataBase.",
                    default = "UniProt")
parser.add_argument("--hmm_db_name", help = "name to be assigned to the hmm database to be created. Its recomended to give a name that \
                    that describes the family or other characteristic of the given sequences")
parser.add_argument("-it", "--input_type", default = "protein", help = "specifies the nature of the sequences in the input file between \
                    'protein', 'nucleic' or 'metagenome'. Defaults to 'protein'")
parser.add_argument("-o", "--output", default = "MPARTY_results", help = "name for the output directory. Defaults to 'MPARTY_results'")
parser.add_argument("--output_type", default = "tsv", help = "choose report table outpt format from 'tsv', 'csv' or 'excel'. Defaults to 'tsv'")
parser.add_argument("-rt", "--report_text", default = False, action = "store_true", help = "decides whether to produce or not a friendly report in \
                    txt format with easy to read information")
parser.add_argument("--hmms_output_type", default = "tsv", help = "chose output type of hmmsearch run from 'out', 'tsv' or 'pfam' format. Defaults to 'tsv'")
parser.add_argument("--validation", default = False, action = "store_true", help = "decides whether to perform models validation and filtration with \
                    the 'leave-one-out' cross validation methods. Call to set to True. Defaults to False")
parser.add_argument("-p", "--produce_inter_tables", default = False, action = "store_true", help = "call if user wants to save intermediate\
                    tables as parseale .csv files (tables from hmmsearch results processing)")
parser.add_argument("--negative_db", help = "path to a user defined negative control database. Default use of human gut microbiome")
parser.add_argument("-s", "--snakefile", help = "user defined snakemake workflow Snakefile. Defaults to '/workflow/Snakefile",
                    default = "workflow/Snakefile")
parser.add_argument("-t", "--threads", type = int, help = "number of threads for Snakemake to use. Defaults to 1",
                    default = 1)
parser.add_argument("--align_method", default = "diamond", help = "chose the alignment method for the initial sequences database expansion, between\
                    'diamond', 'blast' and 'upimapi'. Defaults to 'upimapi'")
parser.add_argument("--aligner", default = "tcoffee", help = "chose the aligner program to perform the multiple sequence alignment for the models\
                    between 'tcoffee' and 'muscle'. Defaults to 'tcoffee'.")
parser.add_argument("-hm", "--hmm_models", type=str, help = f"path to a directory containing HMM models previously created by the user. By default\
                    M-PARTY uses the built-in HMMs from database in 'resources/Data/HMMs/After_tcoffee_UPI/'")
parser.add_argument("--concat_hmm_models", action = "store_false", default = True, help = "call to not concatenate HMM models into a single file. Defaults to True")
parser.add_argument("--unlock", action = "store_true", default = False, help = "could be required after forced workflow termination")
parser.add_argument("-w", "--workflow", default = "annotation", help = 'defines the workflow to follow,\
                    between "annotation", "database_construction" and "both". Latter keyword makes the database construction\
                    first and posterior annotation. Defaults to "annotation"')
parser.add_argument("-c", "--config_file", help = "user defined config file. Only recommended for\
                    advanced users. Defaults to 'config.yaml'. If given, overrides config file construction\
                    from input", default = "config.yaml")
parser.add_argument("--overwrite", action = "store_true", default = False, help = "Call to overwrite inputted files. Defaults to False")
parser.add_argument("--verbose", action = "store_true", default = False, help = "Call so M-PARTY display more messaging")
parser.add_argument("--display_config", default = False, action = "store_true", 
                    help = "declare to output the written config file together with results. Useful in case of debug")
parser.add_argument("-v", "--version", action = "version", version = "M-PARTY {}".format(version))
args = parser.parse_args()
print(vars(args))


strat = "/".join(sys.path[0].split("/")[:-1])
snakefile_path = sys.path[1].replace("\\", "/")+"/workflow/Snakefile"
# config_path = "/".join(sys.path[0].split("\\")[:-1])+"/config/config.yaml"  # for WINDOWS
config_path = "/".join(sys.path[1].split("/"))+"/config/"  # for Linux
hmm_database_path = f'{"/".join(sys.path[1].split("/"))}/resources/Data/HMMs/{args.hmm_db_name}/After_tcoffee_UPI/'
validated_hmm_dir = f'{"/".join(sys.path[1].split("/"))}/resources/Data/HMMs/{args.hmm_db_name}/validated_HMM/'


def read_config_yaml(filename: str) -> tuple:
    config_type = filename.split(".")[-1]
    if config_type == "yaml":
        with open(filename) as stream:
            try:
                config_file = yaml.safe_load(stream)
                stream.close()
            except yaml.YAMLError as exc:
                print(exc)
    elif config_type == "json":
        with open(filename) as stream:
            try:
                config_file == json.load(stream)
                stream.close()
            except json.decoder.JSONDecodeError as exc:
                print(exc)
    else:
        quit("Config file must be in .yaml or .json format! Get an example config file from config/ folder.")
    return config_file, config_type


def parse_fasta(filename: str, remove_excess_ID: bool = True, meta_gen: bool = False) -> list:
    """Given a FASTA file, returns the IDs from all sequences in that file.
    If file not present, program will be quited and TypeError message raised.

    Args:
        filename (str): Name of FASTA file.
        remove_excess_ID (bool, optional): Decide wether to remove the excess part of UniProt IDs. Defaults to True.
        meta_gen(bool, optional): Set to True if input file is from a metagenomic sample. Defaults to False. Metagenomic samples
        usually have mix IDs

    Returns:
        list: A list containing IDs from all sequences
    """
    unip_IDS = []
    # if only validation, input sequences are not needed
    if args.validation == True and args.workflow == "annotation" and args.input == None:
        return unip_IDS
    elif args.workflow == "database_construction" and args.input == None:
        return unip_IDS
    else:
        try:
            with open(filename, "r") as f:
                try:
                    Lines = f.readlines()
                    for line in Lines:
                        if line.startswith(">"):
                            if meta_gen:
                                if not remove_excess_ID:
                                    unip_IDS.append(line.split(" ")[0][1:])
                                else:
                                    try:
                                        identi = re.findall("\|.*\|", line)
                                        identi = re.sub("\|", "", identi[0])
                                        unip_IDS.append(identi)
                                    except:
                                        identi = line.split(" ")[0]
                                        identi = re.sub(">", "", identi)
                                        unip_IDS.append(identi)
                            else:   
                                if not remove_excess_ID:
                                    unip_IDS.append(line.split(" ")[0][1:])
                                else:
                                    identi = re.findall("\|.*\|", line)
                                    identi = re.sub("\|", "", identi[0])
                                    unip_IDS.append(identi)
                except:
                    quit("File must be in FASTA format.")
        except TypeError:
            raise TypeError("Missing input file! Make sure -i option is filled")
        return unip_IDS


def check_results_directory(output: str) -> str:
    """Automatically creats the path where output should appear. Checks if folder already exists or not in the 
    execution path
    Args:
        output (str): Name for the output folder
    Returns:
        str: Path for the output folder
    """
    Path(output).mkdir(exist_ok=True, parents=True)


def write_config(input_file: str, out_dir: str, config_filename: str, with_results: bool = False) -> yaml:
    """Given a input file, output directory, and a name to assign to the new config file, write that same config file
    accordingly to the given arguments

    Args:
        input_file (str): Name for the input FASTA file
        out_dir (str): Name for the output directory where result shall be directed
        config_filename (str): Name to be given to the new config file

    Returns:
        document: Returns a .yaml format config file, with the given arguments though the CLI
    """
    seq_IDS = parse_fasta(input_file, meta_gen = True if args.input_type == "metagenome" else False)
    if args.validation and args.workflow != "database_construction" and args.workflow != "both" and args.input == None:
        out_dir = None
    else:
        check_results_directory(out_dir)
    dict_file = {"seqids": seq_IDS,
                "database": args.database,
                "input_file": None if seq_IDS == [] else args.input.split("/")[-1],
                "input_file_db_const": args.input_seqs_db_const,
                "hmm_database_name": args.hmm_db_name,
                "alignment_method": args.align_method,
                "msa_aligner": args.aligner,
                "input_type": args.input_type,
                "metagenomic": True if args.input_type == "metagenome" else False,
                "validation": True if args.workflow == "database_construction" or 
                args.workflow == "both" or 
                args.validation else False,
                "concat_models": args.concat_hmm_models,
                "output_directory": out_dir,
                "out_table_format": args.output_type,
                "hmmsearch_out_type": args.hmms_output_type,
                "threads": args.threads,
                "workflow": args.workflow,
                "thresholds": ["60-65", "65-70", "70-75", "75-80", "80-85", "85-90"],
                "verbose": args.verbose,
                "overwrite": args.overwrite}
    Path(config_path).mkdir(parents = True, exist_ok = True)
    caminho = config_path + "/" + config_filename
    config_type = config_filename.split(".")[-1]
    if with_results:
        if config_type == "yaml":
            with open(f'{out_dir}/{config_filename}', "w") as file:
                document = yaml.dump(dict_file, file)
                file.close()
        else:
            with open(f'{out_dir}/{config_filename}', "w") as file:
                document = json.dumps(dict_file)
                file.write(document)
                file.close()
    else:
        if config_type == "yaml":
            with open(caminho, "w") as file:
                document = yaml.dump(dict_file, file)
                file.close()
        else:
            with open(caminho, "w") as file:
                document = json.dumps(dict_file)
                file.write(document)
                file.close()
    return document


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


def table_report(dataframe: pd.DataFrame, path: str, type_format: str, db_name: str):
    """Saves a table in a user specified format, with the processed and filtered information from the 
    hmmsearch execution with the HMMs against the query sequences.

    Args:
        dataframe (pd.DataFrame): Dataframe with only the relevant information from hmmsearch execution 
        for all hmm from all threshold ranges.
        path (str): output path.
        type_format (str): Specify the output format.
        db_name(str): Name of the databases name.

    Raises:
        TypeError: Raises TypeError error if user gives an unsupported output format.
    """
    prefix_model = db_name + "_"
    summary_dic = {
        "models": [prefix_model + model for model in get_models_names(dataframe, to_list = True, only_relevant = True)], 
        "querys": get_match_IDS(dataframe, to_list = True, only_relevant = True),
        "bit_scores": get_bit_scores(dataframe, to_list = True, only_relevant = True),
        "e_values": get_e_values(dataframe, to_list = True, only_relevant = True)
        }
    # print(summary_dic)
    indexes = dataframe.index.values.tolist()
    for i in range(len(indexes)):
        summary_dic["models"][i] = f'{indexes[i]}_{summary_dic["models"][i]}'
    df = pd.DataFrame.from_dict(summary_dic)
    table_name = "report_table." + type_format
    if type_format == "tsv":
        df.to_csv(path + table_name, sep = "\t")
    elif type_format == "csv":
        df.to_csv(path + table_name)
    elif type_format == "excel":
        mother_seqs = f'{sys.path[1]}/resources/Data/FASTA/{db_name}/CDHIT/'
        list_IDS_permodel = {}
        for val in summary_dic["models"]:
            thresh = val.split("_")[0]
            model = val.split("_")[-1]
            for folder in os.listdir(mother_seqs):
                if os.path.isdir(os.path.join(mother_seqs, folder)) and folder == thresh:
                    for file in os.listdir(os.path.join(mother_seqs, folder)):
                        if file.endswith(".fasta") and model == file.split(".")[0]:
                            key = thresh + "_" + model
                            if key not in list_IDS_permodel:
                                list_IDS_permodel[key] = parse_fasta(mother_seqs + folder + "/" + file, meta_gen = True if args.input_type == "metagenome" else False)
                            # else:
                            #     list_IDS_permodel[key].append(parse_fasta(mother_seqs + folder + "/" + file, meta_gen = True if args.input_type == "metagenome" else False))
        # print(list_IDS_permodel)
        writer = pd.ExcelWriter(path + "report_table.xlsx", engine = "openpyxl")
        df.to_excel(writer, sheet_name = "Table_Report", index = 0)
        df1 = pd.DataFrame.from_dict(list_IDS_permodel, orient = "index")
        # df1 = df1.transpose()
        df1.to_excel(writer, sheet_name = "Model_Sequences")
        writer.save()  # FutureWarning: save is not part of the public API, usage can give unexpected results and will be removed in a future version
        writer.close()
    else:
        raise TypeError("Specified table format is not available. Read documentation for --output_type.")


def text_report(dataframe: pd.DataFrame, path: str, bit_threshold: float, eval_threshold: float, vali: bool = False):
    """Write the final report as .txt file, with a summary of the results from the annotation 
    performed with hmmsearch. Starts by calculating the number of in-built HMM profiles, and gives an insight of the 
    filtration thresholds.

    Args:
        dataframe (pd.DataFrame): Dataframe with only the relevant information from hmmsearch execution 
        for all hmm from all threshold ranges.
        path (str): output path.
    """
    # number of initial HMM profiles
    number_init_hmms, number_validated_hmms = 0, 0
    for dir in os.listdir(hmm_database_path):
        if os.path.isdir(os.path.join(hmm_database_path, dir)):
            for _ in os.listdir(os.path.join(hmm_database_path, dir)):
                number_init_hmms += 1
    if vali:
        for dir in os.listdir(validated_hmm_dir):
            if os.path.isdir(os.path.join(validated_hmm_dir, dir)):
                for _ in os.listdir(os.path.join(validated_hmm_dir, dir)):
                    number_validated_hmms += 1
    # get the IDs from all hits after quality check
    query_names = get_match_IDS(dataframe, to_list = True, only_relevant = True)
    # get number of hits given for each sequence
    number_hits_perseq = get_number_hits_perseq(query_names)
    # get the unique sequences
    unique_seqs = get_unique_hits(query_names)
    inputed_seqs = parse_fasta(args.input, meta_gen = config["metagenomic"])
    if config["validation"] == True:
        with open(path + "text_report.txt", "w") as f:
            f.write(f"M-PARTY hits report:\n \
                    \nFrom a total number of {number_init_hmms} HMM profiles initially considered, only {len(query_names)} where considered"
                    f"for the final report.\nUser defined validation to true with {args.negative_db} database, from which resulted"
                    f"in {number_validated_hmms}. After annotation, another filtering process was performed considering the values from bit score and E-value from the HMM search run, \n"
                    f"in which the considered bit score threshold was {bit_threshold} and E-value was {eval_threshold}.\n"
                    f"Also, {len(inputed_seqs)} initial query sequences where inputed form {args.input} file, from which {len(unique_seqs)}\n"
                    f"out of these {len(inputed_seqs)} were considered to have a hit against the HMM database.")
            f.close
    else:
        with open(path + "text_report.txt", "w") as f:
            f.write(f"M-PARTY hits report:\n \
                    \nFrom a total number of {number_init_hmms} HMM profiles initially considered, only {len(query_names)} where considered"
                    "for the final report.\nFiltering process was performed considering the values from bit score and E-value from the HMM search run, \n"
                    f"in which the considered bit score threshold was {bit_threshold} and E-value was {eval_threshold}.\n"
                    f"Also, {len(inputed_seqs)} initial query sequences where inputed form {args.input} file, from which {len(unique_seqs)}\n"
                    f"out of these {len(inputed_seqs)} were considered to have a hit against the HMM database.")
            f.close


def get_number_hits_perseq(hit_IDs_list: list) -> dict:
    """Given a list of sequences IDs from the hits against the hmm models from hmmsearch, counts the number of each ID.

    Args:
        hit_IDs_list (list): List of sequence IDs.

    Returns:
        dict: Dictionary containing each ID as key and the respective number of occurrences as value.
    """
    counter = {}
    for i in hit_IDs_list:
        counter[i] = counter.get(i, 0) + 1
    return counter


def get_unique_hits(hit_IDs_list: list) -> list:
    """Given a list of sequence IDs from the hits against the hmm models from hmmsearch, return a new list with only the unique elements.

    Args:
        hit_IDs_list (list): List of sequence IDs.

    Returns:
        list: List with only a single occurrence of each ID.
    """
    unique_IDs_list = []
    for x in hit_IDs_list:
        if x not in unique_IDs_list:
            unique_IDs_list.append(x)
    return unique_IDs_list


def get_aligned_seqs(hit_IDs_list: list, path: str, inputed_seqs: str):
    """Wirtes an ouput Fasta file with the sequences from the input files that had a hit in hmmsearch 
    annotation against the hmm models.

    Args:
        hit_IDs_list (list): list of Uniprot IDs that hit.
        path (str): ouput path.
        inputed_seqs (str): name of the initial input file.
    """
    with open(path + "aligned.fasta", "w") as wf:
        # returns list of IDs from inputed FASTA sequences (entire ID)
        input_IDs = parse_fasta(inputed_seqs, remove_excess_ID = False, meta_gen = config["metagenomic"])
        # print(input_IDs[0:10])
        # print("Sequencias que vieram do input file", input_IDs)
        # time.sleep(5)
        # returns a list the sequences that hit against the models (only one entry)
        unique_IDS = get_unique_hits(hit_IDs_list)
        # print(unique_IDS)
        # print("Sequencias que vieram dos hits com os HMMs", unique_IDS)
        # time.sleep(5)
        with open(inputed_seqs, "r") as rf:
            Lines = rf.readlines()
            for x in unique_IDS:
                if x in input_IDs:
                    # try:
                        iterador = iter(Lines)
                        linha = next(iterador)
                        while linha is not None:
                            if x not in linha:
                                linha = next(iterador, None)
                                continue
                            elif x in linha:
                                wf.write(linha)
                                # print(linha)
                                linha = next(iterador, None)
                                # print(linha)
                                while linha is not None and not linha.startswith(">"):
                                    wf.write(linha)
                                    linha = next(iterador, None)
                                    # print(linha)
                            elif x not in linha and linha.startswith(">"):
                                break
                            linha = next(iterador, None)
                    # except:
                        # quit("File must be in Fasta format.")
                else:
                    continue
        rf.close()
    wf.close()


def generate_output_files(dataframe: pd.DataFrame, hit_IDs_list: list, inputed_seqs: str, bit_threshold: float, eval_threshold: float):
    """Function that initializes the output files creation simultaneously, for now, only two files are generated:
    report and aligned sequences.
    Path will always be the output folder defined by the user when running tool in CLI, so no pat argument is required.

    Args:
        dataframe (pd.DataFrame): Dataframe with only the relevant information from hmmsearch execution.
        hit_IDs_list (list): list of Uniprot IDs that hit.
        inputed_seqs (str): name of the initial input file.
    """
    out_folder = args.output + "/"
    table_report(dataframe, out_folder, args.output_type, args.hmm_db_name)
    if args.report_text:
        if args.validation:
            text_report(dataframe, out_folder, bit_threshold, eval_threshold, vali = True)
        else:
            text_report(dataframe, out_folder, bit_threshold, eval_threshold)
    get_aligned_seqs(hit_IDs_list, out_folder, inputed_seqs)
    if args.display_config:
        write_config(args.input, args.output, args.config_file, with_results = True)

doc = write_config(args.input, args.output, args.config_file)
config, config_format = read_config_yaml(config_path + args.config_file)

if args.hmm_db_name is None:
    raise TypeError("Missing hmm database name! Make sure --hmm_db_name option is filled")
else:
    hmmsearch_results_path = sys.path[1].replace("\\", "/") + "/results/" + args.hmm_db_name + "/HMMsearch_results/"


st = time.time()

# first only runs for if user flags --validation alone without input sequences, will validate the models inside database only
if args.validation and args.workflow != "database_construction" and args.workflow != "both" and args.input == None:

    print("Starting validation procedures...\n")
    time.sleep(2)

    if args.hmm_db_name is None:
        raise TypeError("Missing hmm database name! Make sure --hmm_db_name option is filled")
    else:
        pathing = make_paths_dic(args.hmm_db_name)
    exec_testing(thresholds = config["thresholds"], path_dictionary = pathing, database = args.negative_db)
    to_remove = hmm_filtration(pathing, verbose = config["verbose"])
    remove_fp_models(to_remove, pathing, config["verbose"])
    concat_final_model(pathing)
    time.sleep(2)
    print("M-PARTY has concluded model validation! Will now switch to the newlly created models (in the validated_HMM folder\n")

# runs if input sequences are given
if args.workflow == "annotation" and args.input is not None:

    print("Annotation workflow with hmmsearch started...\n")
    time.sleep(2)

    Path(hmmsearch_results_path).mkdir(parents = True, exist_ok = True)
    if args.validation:
        if not os.path.exists(validated_hmm_dir):
            print("Starting validation procedures...\n")
            time.sleep(2)

            if args.hmm_db_name is None:
                raise TypeError("Missing hmm database name! Make sure --hmm_db_name option is filled")
            else:
                pathing = make_paths_dic(args.hmm_db_name)
            exec_testing(thresholds = config["thresholds"], path_dictionary = pathing ,database = args.negative_db)
            to_remove = hmm_filtration(pathing, verbose = config["verbose"])
            print(to_remove)
            remove_fp_models(to_remove, pathing, config["verbose"])
            concat_final_model(pathing)
            time.sleep(2)
            print("M-PARTY has concluded model validation! Will now switch to the newlly created models (in the validated_HMM folder\n")
        else:
            print("Validated HMM already up, proceding to annotation...\n")
            time.sleep(2)
        if args.concat_hmm_models:
            for hmm_file in file_generator(validated_hmm_dir, full_path = True):
                run_hmmsearch(args.input, hmm_file, 
                            hmmsearch_results_path + "search_" + config["input_file"].split("/")[-1].split(".")[0] +
                            "_" + hmm_file.split("/")[-1].split(".")[0] + "." + args.hmms_output_type,
                            out_type = args.hmms_output_type)
        else:
            p = os.listdir(hmm_database_path)
            for thresh in p:
                path = os.path.join(hmm_database_path, thresh)
                Path(path).mkdir(parents = True, exist_ok = True)
                for hmm_file in file_generator(path, full_path = True):
                    run_hmmsearch(args.input, hmm_file, 
                        path + "search_" + config["input_file"].split("/")[-1].split(".")[0] +
                        "_" + hmm_file.split("/")[-1].split(".")[0] + "." + args.hmms_output_type,
                        out_type = args.hmms_output_type)
                concat_hmmsearch_results(path, hmmsearch_results_path)

    else:
        if args.concat_hmm_models:
            for hmm_file in file_generator(hmm_database_path, full_path = True):
                print(hmm_file)
                run_hmmsearch(args.input, hmm_file, 
                        hmmsearch_results_path + "search_" + config["input_file"].split("/")[-1].split(".")[0] +
                        "_" + hmm_file.split("/")[-1].split(".")[0] + "." + args.hmms_output_type,
                        out_type = args.hmms_output_type)
        else:
            p = os.listdir(hmm_database_path)
            for thresh in p:
                path = os.path.join(hmm_database_path, thresh)
                Path(path).mkdir(parents = True, exist_ok = True)
                for hmm_file in file_generator(path, full_path = True):
                    run_hmmsearch(args.input, hmm_file, 
                        path + "search_" + config["input_file"].split("/")[-1].split(".")[0] +
                        "_" + hmm_file.split("/")[-1].split(".")[0] + "." + args.hmms_output_type,
                        out_type = args.hmms_output_type)
                concat_hmmsearch_results(path, hmmsearch_results_path)

    lista_dataframes = dict.fromkeys(config["thresholds"])
    for file in file_generator(hmmsearch_results_path):
        # if config["input_file"] in file:
            # print(f'File {file} detected \n')
        thresh = file.split("_")[-1].split(".")[0]
        lista_dataframes[thresh] = read_hmmsearch_table(hmmsearch_results_path + file)
    final_df = concat_df_byrow(df_dict = lista_dataframes)
    rel_df = relevant_info_df(final_df)
    # print(rel_df)
    quality_df, bs_thresh, eval_thresh = quality_check(rel_df, give_params = True)
    hited_seqs = get_match_IDS(quality_df, to_list = True, only_relevant = True)
    # print(hited_seqs)
    generate_output_files(quality_df, hited_seqs, args.input, bs_thresh, eval_thresh)

elif args.workflow == "database_construction":
    if args.hmm_db_name is None:
        raise TypeError("Missing hmm database name! Make sure --hmm_db_name option is filled")
    else:
        print("HMM database construction workflow from user input started...\n")
    time.sleep(2)
    if args.input_seqs_db_const is None:
        raise TypeError("Missing input sequences to build HMM database! Make sure --input_seqs_db_const option is filled with a fasta file.")
    time.sleep(2)

    query_DB = build_UPI_query_DB("resources/Data/FASTA/DataBases", config = config, verbose = config["verbose"])

    if config["alignment_method"] == "diamond":
        ### FASTA to DMND
        # diamond_file = build_diamond_DB(args.input_seqs_db_const, "resources/Data/FASTA/", verbose = config["verbose"])
        Path(f'resources/Alignments/{args.hmm_db_name}/BLAST/diamond_output/').mkdir(parents = True, exist_ok = True)
        aligned_TSV = run_UPIMAPI(args.input_seqs_db_const, f'resources/Alignments/{args.hmm_db_name}/BLAST/diamond_output/out.tsv', query_DB, args.threads)
        # aligned_TSV = run_UPIMAPI(query_DB, f'resources/Alignments/{args.hmm_db_name}/BLAST/diamond_output/out.tsv', args.input_seqs_db_const, args.threads)
        handle = UPIMAPI_parser(aligned_TSV)
        dic_enzymes = UPIMAPI_iter_per_sim(handle)
        if config["verbose"]:
            print(f'Saving IDs from the ranges of {config["thresholds"]} percentages of similarity.\n')
        save_as_tsv(dic_enzymes, "resources/Data/Tables/diamond_results_per_sim.tsv")

    elif config["alignment_method"] == "upimapi":
        Path("resources/Data/FASTA/DataBases").mkdir(parents = True, exist_ok = True)
        aligned_TSV = run_UPIMAPI(query_DB, f'resources/Alignments/{args.hmm_db_name}/BLAST/upimapi_results', args.input_seqs_db_const, args.threads)
        handle = UPIMAPI_parser(aligned_TSV)
        dic_enzymes = UPIMAPI_iter_per_sim(handle)
        if config["verbose"]:
            print(f'Saving IDs from the ranges of {config["thresholds"]} percentages of similarity.\n')
        save_as_tsv(dic_enzymes, "resources/Data/Tables/UPIMAPI_results_per_sim.tsv")

    elif config["alignment_method"] == "blast":
        ### FALTA COMANDO PARA CORRER O BLASTP ###
        aligned_TSV = f'resources/Alignments/{args.hmm_db_name}/BLAST/test.tsv'
        handle = BLAST_parser(aligned_TSV)
        dic_enzymes = BLAST_iter_per_sim(handle)
        if config["verbose"]:
            print(f'Saving IDs from the ranges of {config["thresholds"]} percentages of similarity.\n')
        save_as_tsv(dic_enzymes, "resources/Data/Tables/BLAST_results_per_sim.tsv")

    else:
        raise ValueError("--align_method flag only ranges from 'diamond', 'upimapi' or 'blast'. Chose one from the list.")

    from seq_download import get_fasta_sequences

    Path(f'resources/Data/FASTA/{args.hmm_db_name}/UPIMAPI/').mkdir(parents = True, exist_ok = True)
    Path(f'resources/Data/Tables/{args.hmm_db_name}/CDHIT_clusters/').mkdir(parents = True, exist_ok = True)
    for thresh in config["thresholds"]:
        if config["verbose"]:
            print(f'Retrieving sequences from {thresh} range\n')
        get_fasta_sequences("resources/Data/Tables/BLAST_results_per_sim.tsv", f'resources/Data/FASTA/{args.hmm_db_name}/UPIMAPI/{thresh}.fasta')

        ### run CDHIT
        if config["verbose"]:
            print(f'CDHIT run for {thresh} range\n')
        try:
            run_CDHIT(f'resources/Data/FASTA/{args.hmm_db_name}/UPIMAPI/{thresh}.fasta', f'resources/Data/FASTA/{args.hmm_db_name}/UPIMAPI/cd-hit90_after_diamond_{thresh}.fasta')
            handle = cdhit_parser(f'resources/Data/FASTA/{args.hmm_db_name}/UPIMAPI/cd-hit90_after_diamond_{thresh}.fasta.clstr')
            handle2 = counter(handle, tsv_ready = True, remove_duplicates = True)
            save_as_tsv(handle2, f'resources/Data/Tables/{args.hmm_db_name}/CDHIT_clusters/cdhit_clusters_{thresh}_afterUPIMAPI.tsv')

            from CDHIT_seq_download import fasta_retriever_from_cdhit
            if config["verbose"]:
                print(f'Retrieving sequences divided by clusters from CDHIT\n')
            Path(f'resources/Data/FASTA/{args.hmm_db_name}/CDHIT/{thresh}/').mkdir(parents = True, exist_ok = True)
            fasta_retriever_from_cdhit(f'resources/Data/Tables/{args.hmm_db_name}/CDHIT_clusters/cdhit_clusters_{thresh}_afterUPIMAPI.tsv', 
                                        f'resources/Data/FASTA/{args.hmm_db_name}/CDHIT/{thresh}')
        except:
            if config["verbose"]:
                print(f'[WARNING] {thresh} range of similarities not detected.\n')
            time.sleep(2)
            continue

    ### add cluster per threshol to config
    os.remove("config/config.yaml")

    if config_format == "yaml":
        files = get_tsv_files(config)
        threshandclust = threshold2clusters(files)
        print(threshandclust)
        for thresh, cluster in threshandclust.items():
            for c in range(len(cluster)):
                cluster[c] = str(cluster[c])
            config[thresh] = cluster
        newthresh = []
        for thresh in config["thresholds"]:
            if thresh not in threshandclust.keys():
                continue
            else:
                newthresh.append(thresh)
        config["thresholds"] = newthresh
    elif config_format == "json":
        pass

    with open("config/config.yaml", "w") as dump_file:
        yaml.dump(config, dump_file)
        dump_file.close()

    snakemake.main(
        f'-s {args.snakefile} --printshellcmds --cores {config["threads"]} --configfile config/{args.config_file}'
        f'{" --unlock" if args.unlock else ""}')

    files = [f for f in os.listdir('.') if os.path.isfile(f)]
    for file in files:
        if file.endswith(".dnd"):
            delete_inter_files(file)

    print("HMM database created!")
    time.sleep(2)

    if args.validation:
        print("Starting validation procedures...")
        time.sleep(2)
        if args.hmm_db_name is None:
            raise TypeError("Missing hmm database name! Make sure --hmm_db_name option is filled")
        else:
            pathing = make_paths_dic(args.hmm_db_name)
        exec_testing(thresholds = config["thresholds"], path_dictionary = pathing, database = args.negative_db)
        to_remove = hmm_filtration(pathing)
        remove_fp_models(to_remove, pathing)
        concat_final_model(pathing)
        time.sleep(2)
        print("M-PARTY has concluded model validation! Will now switch to the newlly created models (in the validated_HMM folder")

elif args.workflow == "both":
    if args.hmm_db_name is None:
        raise TypeError("Missing hmm database name! Make sure --hmm_db_name option is filled")
    else:
        print("HMM database construction workflow from user input started...")
    time.sleep(2)

    ### UPIMAPI run DIAMOND
    Path("resources/Data/FASTA/DataBases/").mkdir(parents = True, exist_ok = True)
    query_DB = build_UPI_query_DB("resources/Data/FASTA/DataBases", config = config)
    # aligned_TSV = run_UPIMAPI(query_DB, f'resources/Alignments/{args.hmm_db_name}/BLAST/upimapi_results', args.input_seqs_db_const, args.threads)
    aligned_TSV = run_UPIMAPI(query_DB, f'resources/Alignments/{args.hmm_db_name}/BLAST/diamond_results/DIAMOND_results.out', args.input_seqs_db_const, args.threads)
    handle = UPIMAPI_parser(aligned_TSV)
    dic_enzymes = UPIMAPI_iter_per_sim(handle)
    save_as_tsv(dic_enzymes, "resources/Data/Tables/diamond_results_per_sim.tsv")

    from seq_download import get_fasta_sequences

    Path(f'resources/Data/FASTA/{args.hmm_db_name}/UPIMAPI/').mkdir(parents = True, exist_ok = True)
    for thresh in config["thresholds"]:
        get_fasta_sequences("resources/Data/Tables/UPIMAPI_results_per_sim.tsv", f'resources/Data/FASTA/{args.hmm_db_name}/UPIMAPI/{thresh}.fasta')

        ### run CDHIT
        run_CDHIT(f'resources/Data/FASTA/{args.hmm_db_name}/UPIMAPI/{thresh}.fasta', f'resources/Data/FASTA/{args.hmm_db_name}/UPIMAPI/cd-hit90_after_diamond_{thresh}.fasta')
        handle = cdhit_parser(f'resources/Data/FASTA/{args.hmm_db_name}/UPIMAPI/cd-hit90_after_diamond_{thresh}.fasta.clstr')
        handle2 = counter(handle, tsv_ready = True, remove_duplicates = True)
        save_as_tsv(handle2, f'resources/Data/Tables/{args.hmm_db_name}/CDHIT_clusters/cdhit_clusters_{thresh}_afterUPIMAPI.tsv')

        from CDHIT_seq_download import fasta_retriever_from_cdhit
        fasta_retriever_from_cdhit(f'resources/Data/Tables/{args.hmm_db_name}/CDHIT_clusters/cdhit_clusters_{thresh}_afterUPIMAPI.tsv', 
                                    f'resources/Data/FASTA/{args.hmm_db_name}/CDHIT/{thresh}')

    ### add cluster per threshol to config
    os.remove("config/config.yaml")

    if config_format == "yaml":
        files = get_tsv_files(config)
        # print(files)
        threshandclust = threshold2clusters(files)
        # print(threshandclust)
        for thresh, cluster in threshandclust.items():
            for c in range(len(cluster)):
                cluster[c] = str(cluster[c])
            config[thresh] = cluster
        # print(conf_file)
    elif config_format == "json":
        pass

    with open("config/config.yaml", "w") as dump_file:
        yaml.dump(config, dump_file)
        dump_file.close()

    snakemake.main(
        f'-s {args.snakefile} --printshellcmds --cores {config["threads"]} --configfile {args.config_file}'
        f'{" --unlock" if args.unlock else ""}')

    print("HMM database created!")

    files = [f for f in os.listdir('.') if os.path.isfile(f)]
    for file in files:
        if file.endswith(".dnd"):
            delete_inter_files(file)

    time.sleep(2)

    if args.validation:
        print("Starting validation procedures...")
        time.sleep(2)
        if args.hmm_db_name is None:
            raise TypeError("Missing hmm database name! Make sure --hmm_db_name option is filled")
        else:
            pathing = make_paths_dic(args.hmm_db_name)
        exec_testing(thresholds = config["thresholds"], path_dictionary = pathing, database = args.negative_db)
        to_remove = hmm_filtration(pathing)
        remove_fp_models(to_remove, pathing)
        concat_final_model(pathing)
        time.sleep(2)
        print("M-PARTY has concluded model validation! Will now switch to the newlly created models (in the validated_HMM folder")

    print("Annotation workflow with hmmsearch started...")
    time.sleep(2)
    if args.validation:
        for hmm_file in file_generator(validated_hmm_dir, full_path = True):
            run_hmmsearch(args.input, hmm_file, 
                        hmmsearch_results_path + "search_" + config["input_file"].split("/")[-1].split(".")[0] +
                        "_" + hmm_file.split("/")[-1].split(".")[0] + "." + args.hmms_output_type,
                        out_type = args.hmms_output_type)
    else:
        for hmm_file in file_generator(hmm_database_path, full_path = True):
            run_hmmsearch(args.input, hmm_file, 
                        hmmsearch_results_path + "search_" + config["input_file"].split("/")[-1].split(".")[0] +
                        "_" + hmm_file.split("/")[-1].split(".")[0] + "." + args.hmms_output_type,
                        out_type = args.hmms_output_type)
        lista_dataframes = []
        for file in file_generator(hmmsearch_results_path):
            # print(f'File {file} detected \n')
            lista_dataframes.append(read_hmmsearch_table(hmmsearch_results_path + file))
        final_df = concat_df_byrow(list_df = lista_dataframes)
        rel_df = relevant_info_df(final_df)
        # print(rel_df)
        quality_df, bs_thresh, eval_thresh = quality_check(rel_df, give_params = True)
        hited_seqs = get_match_IDS(quality_df, to_list = True, only_relevant = True)
        # print(hited_seqs)
        generate_output_files(quality_df, hited_seqs, args.input, bs_thresh, eval_thresh)

    print("Exiting M-PARTY's program execution...")

elif args.workflow != "annotation" and args.workflow != "database_construction" and args.workflow != "both":
    raise ValueError("-w worflow flag only ranges from 'annotation', 'database_construction' or 'both'. Chose one from the list.")


et = time.time()
elapsed_time = et - st
elapsed_time = elapsed_time * 1000
print(f'Execution time: {elapsed_time:.4f} milliseconds!')
print(f'M-PARTY has stoped running! Results are displayed in the {args.output} folder :)')
