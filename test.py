from pathlib import Path
import sys
import os
print(sys.path)
from workflow.scripts.hmm_process import concat_df_byrow, read_hmmsearch_table, relevant_info_df


def check_results_directory(output: str) -> str:
    """Automatically creats the path where output should appear. Checks if folder already exists or not in the 
    execution path
    Args:
        output (str): Name for the output folder
    Returns:
        str: Path for the output folder
    """
    Path(output).mkdir(exist_ok=True, parents=True)

lista_dataframes = dict.fromkeys(["60-65", "65-70", "70-75", "75-80", "80-85", "85-90"])

hmmsearch_results_path = sys.path[0].replace("\\", "/")+"/results/HMMsearch_results/"

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


for file in file_generator(hmmsearch_results_path):
        # if config["input_file"] in file:
            # print(f'File {file} detected \n')
    thresh = file.split("_")[-1].split(".")[0]
    lista_dataframes[thresh] = read_hmmsearch_table(hmmsearch_results_path + file)
final_df = concat_df_byrow(df_dict = lista_dataframes)

print(relevant_info_df(final_df))

