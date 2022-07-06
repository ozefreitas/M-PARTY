from scripts.tsv_parser import UPIMAPI_parser, UPIMAPI_iter_per_sim
from scripts.tsv_parser import diamond_parser, iter_per_sim, above_60, devide_by_query
from scripts.uniprot_retriever import fasta_retriever, fasta_retriever_from_cdhit
from scripts.docker_run import docker_run_tcoffee, docker_run_hmmbuild, docker_run_hmmsearch
from scripts.CDHIT_parser import cdhit_parser, counter, save_as_tsv
from scripts.CDHIT_parser import get_clusters
from glob import glob
import os
import re
from itertools import product


# df = diamond_parser("C:/Users/jpsfr/OneDrive/Ambiente de Trabalho/TOOL/PDETool/src/PDEFinder/Alignments/Diamond/diamond_out.tsv")

# enzymes = iter_per_sim(df)
# print(enzymes)

# enzyme_search = above_60(df)
# print(enzyme_search)
# fasta_retriever(enzyme_search, filename="Diamond_target_enzymes")
# 
# enzymes = devide_by_query(df)
# print(enzymes)
# fasta_retriever(enzymes, dic=True, filename="Diamond_target_enzymes_byquery", query=True)


### MÉTODO FINAL COM UPIMAPI ###

# upi = UPIMAPI_parser("C:/Users/jpsfr/OneDrive/Ambiente de Trabalho/TOOL/PDETool/workflow/PDEFinder/Alignments/BLAST/results/upimapi_results/UPIMAPI_results.tsv")
# print(upi)
# 
# upi_enzymes = UPIMAPI_iter_per_sim(upi)
# print(upi_enzymes.keys())
# print(upi_enzymes)

# Para o snakemake, os IDS das enzimas que vieram do UPIMAPI e que foram divididas por thresholds, passam a .tsv

# fasta_retriever(upi_enzymes, dic=True)  # download das sequencias


### Fazer cdhit que foi para os ficheiros C:\Users\jpsfr\OneDrive\Ambiente de Trabalho\TOOL\PDETool\workflow\PDEFinder\Data\FASTA\UPIMAPI\cd-hit90_after_diamond_60-65.fasta.clstr

# Para o snakemake, os ficheiros .clstr seráo processados por cluster e transformados em .tsv, por cada threshold

# for perc in range(60, 86, 5):
#     chave = str(perc)+"-"+str(perc+5)
#     sequencias = cdhit_parser("C:/Users/Ze/Desktop/Mestrado/3ºSemestre/TESE/PDETool/workflow/Data/FASTA/UPIMAPI/cd-hit90_after_diamond_" + chave + ".fasta.clstr")
#     numb_seq = counter(sequencias, remove_single=True, tsv_ready=True)
#     save_as_tsv(numb_seq, "C:/Users/Ze/Desktop/Mestrado/3ºSemestre/TESE/PDETool/workflow/Data/Tables/cdhit_clusters_" + chave + "_afterUPIMAPI.tsv")

### Processamento dos cd-hit depois do upimapi e criacao dos ficheiros fasta com cada cluster

# for perc in range(60, 86, 5):
#     chave = str(perc)+"-"+str(perc+5)
#     sequencias = cdhit_parser("C:/Users/jpsfr/OneDrive/Ambiente de Trabalho/TOOL/PDETool/src/PDEFinder/Data/FASTA/UPIMAPI/cd-hit90_after_diamond_" + chave + ".fasta.clstr")
#     numb_seq = counter(sequencias, remove_single=True)
#     print(numb_seq)
#     fasta_retriever_from_cdhit(numb_seq, chave)


thresholds =  ["60-65", "65-70", "70-75", "75-80", "80-85", "85-90"]

files = {threshold: glob(f"workflow/Data/Tables/cdhit_clusters_{threshold}_afterUPIMAPI.tsv") for threshold in thresholds}
threshold2clusters = {}
for thresh, path in files.items():
	threshold2clusters[thresh] = get_clusters(path[0])
print(files)

# fazer uma lista de listas com todos os clusters, por ordem de threshold
big_list_clusters = [v for k, v in threshold2clusters.items()]
max_clusters = max([max(x) for x in big_list_clusters])
all_clusters = [str(i) for i in range(0, max_clusters+1)]
# print(all_clusters)

# função vai fazer todas as combinações entre thresholds e clusters correspondentes
def util(lista_thresholds, lista_de_listas_clusters):
    autorized_combs = []
    for threshold in range(len(lista_thresholds)):
        for cluster in lista_de_listas_clusters[threshold]:
            combinacao = (lista_thresholds[threshold], str(cluster))
            autorized_combs.append(combinacao)
    autorized_combs_frozen = {frozenset(t) for t in autorized_combs}
    return autorized_combs_frozen

# função que vai buscar os clusters correspondentes a cada threshold
def match_threshold_W_cluster(combinador, desired_combs) -> tuple:
    def match_threshold_W_cluster(*args, **kwargs):
        for combo in combinador(*args, **kwargs):
            if frozenset(combo) in desired_combs:
                yield combo
    return match_threshold_W_cluster


# desired = util(thresholds, big_list_clusters)
# print(desired)
# produto_cartesiano = product(thresholds, all_clusters)
# sucess = 0
# for x in produto_cartesiano:
#     if frozenset(x) in desired:
#        print(x)


files = {threshold: glob(f"workflow/Data/FASTA/CDHIT/{threshold}/*.fasta") for threshold in thresholds}
print(files)
threshold2clusters = {k : [v.split("/")[-1].split("\\")[-1].split('.f')[0] for v in values] for k, values in files.items()}
print(threshold2clusters)

# filtered_product = match_threshold_W_cluster(product, desired)

### Correr t_coffee para todos os ficheiros gerados em cima

# docker_run_tcoffee("C:/Users/jpsfr/OneDrive/Ambiente de Trabalho/TOOL/PDETool/src/PDEFinder/Data/FASTA/CDHIT:/data", "60-65", "1.fasta", "clustalw_aln", "1.clustal_aln")
# docker_run_tcoffee("C:/Users/jpsfr/OneDrive/Ambiente de Trabalho/TOOL/PDETool/src/PDEFinder:/data", "Data/FASTA/CDHIT/60-65", "1.fasta", "clustalw_aln", "Alignments/MultipleSequencesAlign/T_Coffee_UPI/1.clustal_aln")

# cdhit_path = "C:/Users/jpsfr/OneDrive/Ambiente de Trabalho/TOOL/PDETool/src/PDEFinder/Data/FASTA/CDHIT"
# dest_path = "C:/Users/jpsfr/OneDrive/Ambiente de Trabalho/TOOL/PDETool/src/PDEFinder/Alignments/MultipleSequencesAlign/T_Coffee_UPI"
# 
# for fold in os.listdir(cdhit_path):
#     newpath = os.path.join(dest_path, fold)
#     if not os.path.exists(newpath):
#         os.mkdir(newpath)
#     for file in os.listdir(os.path.join(cdhit_path, fold)):
#         print("Executing T-Coffee for file", file, "in the", fold, "thresold.", "\n")
#         docker_run_tcoffee("C:/Users/jpsfr/OneDrive/Ambiente de Trabalho/TOOL/PDETool/src/PDEFinder:/data", "Data/FASTA/CDHIT/" + fold, file, "clustalw_aln", "Alignments/MultipleSequencesAlign/T_Coffee_UPI/" + fold + "/" + file + ".clustal_aln")


### Correr hmmbuild para todos os ficheiros de alinhamento gerados em cima

# tcoffee_path = "C:/Users/jpsfr/OneDrive/Ambiente de Trabalho/TOOL/PDETool/src/PDEFinder/Alignments/MultipleSequencesAlign/T_Coffee_UPI"
# dest_path = "C:/Users/jpsfr/OneDrive/Ambiente de Trabalho/TOOL/PDETool/src/PDEFinder/Data/HMMs/After_tcoffee_UPI"
# 
# for fold in os.listdir(tcoffee_path):
#     newpath = os.path.join(dest_path, fold)
#     if not os.path.exists(newpath):
#         os.mkdir(newpath)
#     for file in os.listdir(os.path.join(tcoffee_path, fold)):
#         if file.endswith(".clustalw_aln"):
#             print("Executing hmmbuild for file", file, "in the", fold, "thresold.", "\n")
#             outname = ""
#             for i in file:
#                 if i != ".":
#                     outname += i
#                 else:
#                     outname += ".hmm"
#                     break
#             print("Output name:", outname)
#             docker_run_hmmbuild("C:/Users/jpsfr/OneDrive/Ambiente de Trabalho/TOOL/PDETool/src/PDEFinder/:/data", "Alignments/MultipleSequencesAlign/T_Coffee_UPI/" + fold + "/" + file, "Data/HMMs/After_tcoffee_UPI/" + fold + "/" + outname)


### Correr hmmsearch para todos os modelos contra a UniProt

# hmmpath = "C:/Users/jpsfr/OneDrive/Ambiente de Trabalho/TOOL/PDETool/src/PDEFinder/Data/HMMs/After_tcoffee_UPI"
# dest_path = "C:/Users/jpsfr/OneDrive/Ambiente de Trabalho/TOOL/PDETool/src/PDEFinder/Data/HMMsearch_results/After_UPI"
# 
# for fold in os.listdir(hmmpath):
#     newpath = os.path.join(dest_path, fold)
#     if not os.path.exists(newpath):
#         os.mkdir(newpath)
#     for file in os.listdir(os.path.join(hmmpath, fold)):
#         print("Executing hmmbuild for file", file, "in the", fold, "thresold.", "\n")
#         outname = ""
#         for i in file:
#             if i != ".":
#                 outname += i
#             else:
#                 outname += ".out"
#                 break
#         print("Output name:", outname)
#         docker_run_hmmsearch("C:/Users/jpsfr/OneDrive/Ambiente de Trabalho/TOOL/PDETool/src/PDEFinder:/data", "Data/HMMs/After_tcoffee_UPI/" + fold + file, 
#                             "C:/Users/jpsfr/OneDrive/Ambiente de Trabalho/TOOL/PDETool/workflow/PDEFinder/Data/FASTA/DataBases/familiesDB.fasta", "Data/HMMsearch_results/After_UPI/" + fold + "/" + outname)
