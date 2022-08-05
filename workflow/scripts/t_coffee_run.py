from docker_run import docker_run_tcoffee

docker_run_tcoffee("/mnt/c/Users/jpsfr/OneDrive/Ambiente de Trabalho/TOOL/PDETool/:/data/", snakemake.input[0], "clustal_aln", snakemake.output[0])