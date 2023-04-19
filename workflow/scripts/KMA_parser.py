from command_run import run_command
import re
import os
import pandas as pd


def run_KMA(input_DB: str, output_DB: str, meta_input: str, meta_out: str, threads: int, paired_end: bool = False, second_input: str = None):
    run_command(f'kma`index`-i`{input_DB}`-o`{output_DB}')
    run_command(f'kma`-i`{meta_input}`{meta_out}`-t_db`{output_DB}`-t`{threads}`-1t1`-mem_mode`-ef')
    if paired_end and second_input == None:
        raise ValueError("When paired end is flaged, a second input file is mandatory")
    elif paired_end and second_input != None:
        run_command(f'kma`-i`{meta_input}`{second_input}`{meta_out}`-t_db`{output_DB}`-t`{threads}`-1t1`-mem_mode`-ef')
    return meta_out

def KMA_parser():
    ficheiro_final = open("C:/Users/diogo/Desktop/ARG_CR/Thesis/temp/kma_" + datab + "_list.txt", "w")
    for id in lista_ids():
        DIR_sample = DIR + datab + "/" + id
        try:
            with open(DIR_sample + "/" + id + "_kma.res", 'r') as f:
                df = pd.read_csv(f, sep='\t')
                df["Hit_len"] = (df["Template_length"] * df["Template_Coverage"]) / 100
                df = df[(df['Hit_len'] >= hit_len)]
        except:
            continue
        delimiter = '\n'
        ficheiro_final.write("$" + id + delimiter)
        for name in df.iloc[:, 0]:
            ficheiro_final.write('>' + name + delimiter)
    ficheiro_final.close()
    pass

