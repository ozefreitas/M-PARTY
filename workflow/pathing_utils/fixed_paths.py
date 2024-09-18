import sys

class PathManager:
    snakefile_path = None
    config_path = None
    hmm_database_path = None
    validated_hmm_dir = None

def declare_fixed_paths(args: dict):
    PathManager.snakefile_path = sys.path[0].replace("\\", "/") + "/workflow/Snakefile"
    PathManager.config_path = sys.path[0] + "/config/"  # for Linux
    PathManager.hmm_database_path = f'{sys.path[0]}/resources/Data/HMMs/{args.hmm_db_name}/'
    PathManager.validated_hmm_dir = f'{"/".join(sys.path[0].split("/"))}/resources/Data/HMMs/{args.hmm_db_name}/validated_HMM/'