import pandas as pd
import urllib.request


def fasta_retriever_from_cdhit(tsv_file: str, out_file: str):
    """Given a .tsv file for each similarity threshold with UniProt ID's as data, first column as number of the cluster, 
    write a fasta file with all fasta sequences of the corresponding ID's, by cluster.

    Args:
        tsv_file (string): String containing the name of the .tsv file to be processed.
        out_file (string): Path to the resulting fasta file.
    """
    df = pd.read_csv(tsv_file, sep="\t", index_col=0)
    # numb_cluster = get_number_clusters(tsv_file)
    threshold = out_file.split("/")[-2]
    cluster = out_file.split("/")[-1].split(".f")[0]
    # print(tsv_file)
    # print(out_file)
    # print("cluster", cluster)
    # print("thresh", threshold)
    for index, content in df.iterrows():
        # print(index, type(index))
        if index == int(cluster):
            # print("sucesso")
            # abre o ficheiro no modo write
            #file = open("c:/Users/jpsfr/OneDrive/Ambiente de Trabalho/TOOL/PDETool/workflow/Data/FASTA/CDHIT" + threshold + "/" + str(index) + ".fasta", mode="w")
            file = open(out_file, mode="w")
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
            break

fasta_retriever_from_cdhit(snakemake.input[0], snakemake.output[0])
