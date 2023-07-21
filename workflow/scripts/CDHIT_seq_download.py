import pandas as pd
import urllib.request
from mparty_util import get_clusters, get_number_clusters


def fasta_retriever_from_cdhit(tsv_file: str, out_dir: str):
    """Given a .tsv file for each similarity threshold with UniProt ID's as data, first column as number of the cluster, 
    write a fasta file with all fasta sequences of the corresponding ID's, by cluster.

    Args:
        tsv_file (str): String containing the name of the .tsv file to be processed.
        out_dir (str): Path to the output directory.
    """
    df = pd.read_csv(tsv_file, sep="\t", index_col=0)
    clusters = get_clusters(tsv_file)
    for index, content in df.iterrows():
        # print(index, type(index))
        # abre o ficheiro no modo write
        out_file = out_dir + "/" + str(index) + ".fasta"
        file = open(out_file, mode = "w")
        for seq in list(content):
            try:
                # faz o download da sequencia em formato fasta
                data = urllib.request.urlopen("http://www.uniprot.org/uniprot/" + seq + ".fasta")
            except:
                continue
            dados = data.read()
            # urllib.request.urlopen retorna os dados em forma de bytes, tem que ser convertido em string
            encoding = 'utf-8'
            fasta = dados.decode(encoding)
            fasta = fasta.split("\n")
            # para cada elemento da lista gerado pelo split
            for line in fasta:
                # vai adicionar uma linha ao ficheiro, juntamente com o \n no final, para fazer uma nova linha
                file.write(line)
                file.write("\n")
        file.close()

# fasta_retriever_from_cdhit(snakemake.input[0], snakemake.output[0])
