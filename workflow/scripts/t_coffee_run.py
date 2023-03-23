from command_run import docker_run_tcoffee
import sys

path = sys.path[0]

docker_run_tcoffee(f'{path}:/data/', snakemake.input[0], "clustal_aln", snakemake.output[0])