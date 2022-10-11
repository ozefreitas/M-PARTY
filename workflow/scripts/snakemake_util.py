from glob import glob
from itertools import product
import pandas as pd


def get_clusters(tsv_file):
    df = pd.read_csv(tsv_file, sep="\t", index_col=0)
    return list(df.index.values)


def get_number_clusters(tsv_file):
    df = pd.read_csv(tsv_file, sep="\t", index_col=0)
    return len(df.index)


# vai buscar aos .tsv criados antes, e não fica dependente dos fasta que vão ser criados
def get_tsv_files(config_file):
	files = {threshold: glob(f"workflow/Data/Tables/cdhit_clusters_{threshold}_afterUPIMAPI.tsv") for threshold in config_file["thresholds"]}
	return files


def threshold2clusters(file_list):
	threshold2clusters = {}
	for thresh, path in file_list.items():
		threshold2clusters[thresh] = get_clusters(path[0])


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


def download_uniprot():
	pass

