import yaml
import json
from snakemake_util import threshold2clusters, get_tsv_files
import os


def read_config_yaml(filename: str) -> tuple:
    config_type = filename.split(".")[-1]
    if config_type == "yaml":
        with open(filename) as stream:
            try:
                config_file = yaml.safe_load(stream)
                stream.close()
            except yaml.YAMLError as exc:
                print(exc)
    elif config_type == "json":
        with open(filename) as stream:
            try:
                config_file == json.load(stream)
                stream.close()
            except json.decoder.JSONDecodeError as exc:
                print(exc)
    else:
        quit("Config file must be in .yaml or .json format! Get an example config file from config/ folder.")
    return config_file, config_type

if os.path.exists(snakemake.output[0]):
    os.remove(snakemake.output[0])

conf_file, file_type = read_config_yaml(snakemake.input[0])
# print(conf_file)

os.remove(snakemake.input[0])

if file_type == "yaml":
    files = get_tsv_files(conf_file)
    # print(files)
    threshandclust = threshold2clusters(files)
    # print(threshandclust)
    for thresh, cluster in threshandclust.items():
        for c in range(len(cluster)):
            cluster[c] = str(cluster[c])
        conf_file[thresh] = cluster
    # print(conf_file)
elif file_type == "json":
    pass


with open(snakemake.output[0], "w") as dump_file:
    yaml.dump(conf_file, dump_file)
    dump_file.close()


output = snakemake.output[0]
new_output = output.split("_")[1]
new_output2 = "config/" + new_output


with open(new_output2, "w") as dump_file:
    yaml.dump(conf_file, dump_file)
    dump_file.close()
