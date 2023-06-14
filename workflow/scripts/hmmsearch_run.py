from command_run import run_command

def run_hmmsearch(sequences_file: str, hmm_file: str, output_file: str, verbose: bool = True, eval: float = 10.0, out_type = "out"):
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