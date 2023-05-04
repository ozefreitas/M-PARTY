from command_run import run_command
import re
import os
import pandas as pd


def run_KMA(input_DB: str, output_DB: str, meta_input: str, meta_out: str, threads: int, paired_end: bool = False, second_input: str = None):
    """Runs KMA with in the command line and returns the output file

    Args:
        input_DB (str): _description_
        output_DB (str): _description_
        meta_input (str): _description_
        meta_out (str): _description_
        threads (int): _description_
        paired_end (bool, optional): _description_. Defaults to False.
        second_input (str, optional): _description_. Defaults to None.

    Raises:
        ValueError: _description_

    Returns:
        _type_: _description_
    """
    run_command(f'kma`index`-i`{input_DB}`-o`{output_DB}')
    run_command(f'kma`-i`{meta_input}`-o`{meta_out}`-t_db`{output_DB}`-t`{threads}`-1t1`-mem_mode`-ef')
    if paired_end and second_input == None:
        raise ValueError("When paired end is flaged, a second input file is mandatory")
    elif paired_end and second_input != None:
        run_command(f'kma`-i`{meta_input}`{second_input}`{meta_out}`-t_db`{output_DB}`-t`{threads}`-1t1`-mem_mode`-ef')
    return meta_out

def kma_parser(input_file: str):
    """_summary_

    Args:
        input_file (str): _description_

    Returns:
        _type_: _description_
    """
    df = pd.read_csv(input_file, sep = "\t", decimal = ",")
    df1 = df.loc[(df["Template_Identity"] >= 80) & 
                 (df["Template_Coverage"] >= 60)]
    return df1[["#Template", "Template_Identity", "Template_Coverage", "q_value", "p_value"]]

def get_hit_sequences(dataframe: pd.DataFrame, to_list: bool = False) -> list:
    if to_list:
        return dataframe["#Template"].tolist()
    else:
        return dataframe["#Template"]


# df = kma_parser("/mnt/c/Users/Ze/Desktop/BOLSA/resultados_KMA/SRR3933361_nuc_fastq_kma.res")
# print(df)