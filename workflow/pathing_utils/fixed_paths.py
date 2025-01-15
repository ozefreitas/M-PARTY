import sys
from pathlib import Path


class PathManager:
    system_path = Path(sys.path[0])
    snakefile_path = None
    config_path = None
    hmm_database_path = None
    tcoffee_path = None
    clusters_path = None
    validated_hmm_dir = None
    fasta_type_dir = None

def declare_fixed_paths(args: dict):
    PathManager.snakefile_path = PathManager.system_path / "workflow" / "Snakefile"
    PathManager.config_path = PathManager.system_path / "config"
    PathManager.fasta_type_dir = PathManager.system_path / "resources "/ "Data" / "FASTA" / args.hmm_db_name
    PathManager.hmm_database_path = PathManager.system_path / "resources" / "Data" /" HMMs"/ args.hmm_db_name
    PathManager.tcoffee_path = PathManager.system_path / "resources" / "Alignments" / args.hmm_db_name / "MultipleSequencesAlign" /"T_Coffee"
    PathManager.clusters_path = PathManager.system_path / "resources" / "Data" / "FASTA" / args.hmm_db_name / "CDHIT" / "clusters"
    PathManager.validated_hmm_dir = PathManager.system_path / "resources" / "Data" / "HMMs" / args.hmm_db_name / "validated_HMM"