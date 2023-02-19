import pandas as pd
import re
import urllib.request


def get_fasta_sequences(tsv_file: str, out_folder: str, seq_proc: bool = False):
    """Given a .tsv file with UniProt ID's as data, first column as thresholds, write a fasta file with all fasta sequences of the corresponding ID's,
        by similarity thresholds.
    Args:
        tsv_file (string): String containing the name of the .tsv file to be processed.
        out_folder (string): String regarding the path to the folder where files will be saved.
        seq_proc (bool, optional): Decides if is necessary to process UniProt ID's. Defaults to False.
    """
    # abrir tsv como DataFrame
    df = pd.read_csv(tsv_file, sep="\t", index_col=0)
    threshold = out_folder.split("/")[-1].split(".")[0]
    # print("teste de threshold:", threshold)
    # iterar pelas linhas
    for index, content in df.iterrows():
        if str(index) == threshold:
            file = open(file = "/".join(out_folder.split("/")[:-1]) + "/" + str(index) + ".fasta", mode="w")
            # file = "/".join(out_folder.split("/")[:-1]) + "/" + index + ".fasta"
            for uni_id in list(content):
            # muda a seq para o codigo que o uniprot aceite como ID
                if seq_proc:
                    code = re.findall("\|.*\|", uni_id)
                    clean = re.sub("\|", "", code[0])
                else:
                    clean = uni_id
                try:
                   # faz o download da sequencia em formato fasta
                    data = urllib.request.urlopen("http://www.uniprot.org/uniprot/" + clean + ".fasta")
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

# get_fasta_sequences(snakemake.input[0], snakemake.output[0])