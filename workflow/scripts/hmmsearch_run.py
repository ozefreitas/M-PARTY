from workflow.scripts.command_run import run_command

def run_hmmsearch(sequences_file: str, hmm_file: str, output_file: str, verbose: bool = True, eval: float = 10.0, out_type = "out"):
    """Simplified run function for hmmsearch.

    Args:
        sequences_file (str): file with the sequences to be queried against the models.
        hmm_file (str): file containing the HMMs.
        output_file (str): file where the output will be written.
        verbose (bool, optional): decides to print more information. Defaults to True.
        eval (float, optional): e-value. Defaults to 10.0.
        out_type (str, optional): hmmsearch suports three types of outputs: "out", "tsv" and "pfam". Defaults to "out". This parameter 
        can be changed through user input.
    """
    if out_type == "out":
        run_command(f"hmmsearch`-E`{eval}`{hmm_file}`{sequences_file}`>`{output_file}", sep = "`")
        if not verbose:
            run_command(f"hmmsearch`--noali`-E`{eval}`{hmm_file}`{sequences_file}`>`{output_file}", sep = "`")
    elif out_type == "tsv":
        run_command(f"hmmsearch`-E`{eval}`--tblout`{output_file}`{hmm_file}`{sequences_file}", sep = "`")
        if not verbose:
            run_command(f"hmmsearch`--noali`-E`{eval}`--tblout`{output_file}`{hmm_file}`{sequences_file}", sep = "`")
    elif out_type == "pfam":
        run_command(f"hmmsearch`-E`{eval}`--pfamtblout`{output_file}`{hmm_file}`{sequences_file}", sep = "`")
        if not verbose:
            run_command(f"hmmsearch`--noali`-E`{eval}`--pfamtblout`{output_file}`{hmm_file}`{sequences_file}", sep = "`")