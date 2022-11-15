from glob import glob
from itertools import product
from clint.textui import progress
import pandas as pd
import requests
import shutil
from docker_run import run_command


def get_clusters(tsv_file):
    df = pd.read_csv(tsv_file, sep="\t", index_col=0)
    return list(df.index.values)


def get_number_clusters(tsv_file):
    df = pd.read_csv(tsv_file, sep="\t", index_col=0)
    return len(df.index)


# vai buscar aos .tsv criados antes, e não fica dependente dos fasta que vão ser criados
def get_tsv_files(config_file) -> dict:
	files = {threshold: glob(f"resources/Data/Tables/PE/CDHIT_clusters/cdhit_clusters_{threshold}_afterUPIMAPI.tsv") for threshold in config_file["thresholds"]}
	return files


def threshold2clusters(file_dic: dict) -> dict:
	threshold2clusters = {}
	for thresh, path in file_dic.items():
		threshold2clusters[thresh] = get_clusters(path[0])
	return threshold2clusters


def clusters_in_list(dic: dict) -> list:
	lista_clusters = [v for k, v in dic.items()]
	return lista_clusters


def get_all_clusters(config_file):
# fazer uma lista de listas com todos os clusters, por ordem de threshold
	big_list_clusters = [v for k, v in threshold2clusters.items()]
	max_clusters = max([max(x) for x in big_list_clusters])
	all_clusters = [str(i) for i in range(0, max_clusters+1)]
	return big_list_clusters, all_clusters, max_clusters


# função vai fazer todas as combinações entre thresholds e clusters correspondentes
def util(lista_thresholds, lista_de_listas_clusters):
	autorized_combs = []
	for threshold in range(len(lista_thresholds)):
		for cluster in lista_de_listas_clusters[threshold]:
			combinacao = (lista_thresholds[threshold], str(cluster))
			autorized_combs.append(combinacao)
	autorized_combs_frozen = {frozenset(t) for t in autorized_combs}
	return autorized_combs_frozen


# função que vai fazer o produto entre todos, e so devolve os desejados
def match_threshold_W_cluster(combinador, desired_combs) -> tuple:
    def match_threshold_W_cluster(*args, **kwargs):
        for combo in combinador(*args, **kwargs):
            if frozenset(combo) in desired_combs:
                yield combo
    return match_threshold_W_cluster


# desired = util(config["thresholds"], big_list_clusters)
# inicializar função de combinação
# filtered_product = match_threshold_W_cluster(product, desired)


# def cat_hmms_input(wildcard, config_file):
# 	list_clusters = get_all_clusters(config_file)[0]
# 	return ["workflow/Data/HMMs/After_tcoffee_UPI/{threshold}/{cluster}.hmm".format(threshold=config_file["thresholds"][x], 
# 			cluster=list_clusters[x][y]) for x in range(len(config_file["thresholds"])) for y in range(len(list_clusters[x]))]


def cat_hmms_input(wildcards):
	return expand("resources/Data/HMMs/After_tcoffee_UPI/{threshold}/{cluster}.hmm", threshold=wildcards, cluster=threshold2clusters[wildcards])


def get_target_db(config):
	return config["hmm_database_name"]


def get_output_dir(path, config, hmm = False):
	c = path.split("_")
	if hmm:
		ind = c.index("HMMs")
	else:
		try:
			ind = c.index("FASTA")
		except:
			ind = c.index("Tables")
	c.insert(ind + 1, config["hmm_database_name"])
	return "/".join(c)
 

def get_UPI_queryDB(config):
	return config["database"]


def download_with_progress_bar(url: str, database_folder: str):
	r = requests.get(url, stream=True)
	path = f'{database_folder}/{url.split("/")[-1]}'
	with open(path, "wb") as wf:
		total_length = int(r.headers.get('content-length'))
		for chunk in progress.bar(r.iter_content(chunk_size=1024), expected_size=(total_length/1024) + 1): 
			if chunk:
				wf.write(chunk)
				wf.flush()
	wf.close()
	

def download_uniprot(database_folder: str):
	for url in [
	"https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz", 
	"https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_trembl.fasta.gz",
	"https://ftp.uniprot.org/pub/databases/uniprot/relnotes.txt"]:
		download_with_progress_bar(url, database_folder)
	run_command(f'zcat {database_folder}/{url[0].split("/")[-1]} {database_folder}/{url[1].split("/")[-1]} > {database_folder}/uniprot.fasta')


def build_UPI_query_DB(database_folder, config = None):
	# database = "uniprot"
	database = get_UPI_queryDB(config)
	if database.lower() == "uniprot":
		download_uniprot(database_folder)
	elif database.lower() == "swissprot":
		download_with_progress_bar("https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz", database_folder)
		run_command(f'gunzip -v {database_folder}/uniprot_sprot.fasta.gz')
	elif database.split(".")[-1] == "fasta":
		shutil.move(database, database_folder)
	else:
		raise TypeError("--database given parameter is not accepted. Chose between 'uniprot', 'swissprot' or a path to a FASTA file of protein sequences.")
