# M-PARTY - Mining Protein dAtabases foR Targueted enzYmes
#### Version 0.2.2


<br>
M-PARTY is a free to use, open source user friendly CLI (early release) implemented workflow and database for the detection of plastic degrading enzymes in metagenomic samples, through structural annotation using Hidden Markov Models.
<p>
<br>

## Index

1. [Introduction](https://github.com/ozefreitas/M-PARTY#introduction)
2. [Installation](https://github.com/ozefreitas/M-PARTY#installation)
3. [Usage](https://github.com/ozefreitas/M-PARTY#usage)
4. [Output](https://github.com/ozefreitas/M-PARTY#output)
5. [Additional arguments](https://github.com/ozefreitas/M-PARTY#additional-arguments)

<br>

## Introduction 

M-PARTY is a free to use, open source user friendly CLI implemented workflow and database for the detection of plastic degrading enzymes in metagenomic samples, through structural annotation using Hidden Markov Models, that allows the user to freely interacte with the tool in-built databases and backbone. <p>
M-PARTY compreends a extensive HMM database, built with state of the art checked enzymatic sequences able to degrade plastic polymers, which is used to carry the structural annotation of given sequences. <p>
First version of M-PARTY is only available for mining PE (polyethylene), as latter version will compreend a more vast list of plastics to analyse. <p>
Also, M-PARTY is meant to analyse metagenomic sequences, but version 0.1.0 only accepts single FASTA aminoacidic sequences. Basic steps of M-PARTY annotation workflow in its frist stages are: 

1. The acceptence of any number of protein sequences in a single FASTA file as query;
2. Execution of hmmsearch from the [HMMER package](https://www.hmmer.org/), using the pre-built HMMs from previously knowns sequences able to have some kind of PE deterioration levels as database; 
3. A quality benchmark to determine good and bad hits from the queries against models;
4. Three output files, consisting in a FASTA file with the protein sequences returned as a hit from the search, a report in text format (if requested by the user) with simply puted information about the inputed and already built (HMMs) data, run and processing parameters and conclusions, and an easy to read report table in xlsx format, with all the important data about the annotation results, in particular:
    - Sequence IDs
    - HMM IDs (Degraded plastic + number)
    - Bit scores
    - E-values

<br>

## Installation

M-PARTY is, avaliable for Linux platforms though GitHub repository clonning, using the following line in a git bash terminal inside the desired (empty) folder:

```
cd path/to/desired/dir
git clone https://github.com/ozefreitas/M-PARTY.git 
```

I highly recommed users to create an appropriate conda environment with the required dependencies so M-PARTY executes smoothly, with:

```
cd workflow/envs/ 
conda env create -n <name of env> -f mparty.yaml 
conda activate <name of env> 
cd ../..
```

Clonning though GitHub is only recommended in last case scenario, as this as deprecated in detriment of bioconda distribution aplication.

M-PARTY is available as a conda package from bioconda. Due to tool recent name change, the package still remains with the old name.
Simply open an Anaconda prompt or a command line interface with Anaconda or Miniconda distributions installed and:

```
conda install -c conda-forge -c bioconda m-party
```

and you will be good to go. <p>

If somethig goes wrong, I sugest you to first create a conda environment with:

```
conda create -n <name of env> -c conda-forge -c bioconda m-party
```
due to possible compatibility issues that may occur.

<br>

## Usage

The main and most basic use for M-PARTY is:<p>

```
m-party.py -i path/to/input_file -o output_folder -rt --output_type excel --hmm_db_name PE
```

where the **`-i` input file** must be in FASTA format and contain only (for the time being) aminoacidic sequences, otherwise, program will exit. **`-o` output folder** can be a pre-existing folder or any name for a folder that will be created anyways. The **`-rt`** option flag instructs the tool to include in the output the report in text format, for an easier interpretation of the annotation results and conclusion taking. Also, **`--output_type`** is recommended to be set to "excel" on these earlier versions, as other output format for the table report will be incrementally coded.

<br>

The HMM database M-PARTY has from the start are not validated, has some can be false positive and give deceiving results. If your goal is to only validate a set a HMMs, than run:

```
m-party.py --validation
```

<br>

Now, if you want you can instantly run the annotation workflow from a set of proteins of your liking, and so performing the validation beforehand, with:

```
m-party -i path/to/input_file -o path/to/output_folder -rt --output_type excel --hmm_db_name PE --validation --display_config
```

## Output

M-PARTY will result in three distinct outputs: **report table**, **text report** and **aligned**. In earlier versions, **report table** is only available in *excel* format, although later will also be for *tsv* and *csv*. **Text report** is a user friendly easy to understand summary of the annotation run performed by M-PARTY, and embrace a series of useful information for the user, depending on the given arguments. For last, **aligned** is a FASTA file with all the sequences that had a match in one or more models (this will be refined as model benchmarking and validation are introduced into M-PARTY). 

<br>

## Aditional arguments

M-PARTY is currently in development stage, so only the **"annotation"** workflow is available, and so, also some other features and parameters still have no use and impact in this tool execution, so I highly recommend you to follow the steps in the [usage](https://github.com/ozefreitas/M-PARTY#usage) section.

```
usage: m-party.py [-h] [-i INPUT]
                    [--input_seqs_db_const INPUT_SEQS_DB_CONST]
                    [-ip INPUT_TYPE] [-o OUTPUT] [--output_type OUTPUT_TYPE]
                    [-rt] [--hmms_output_type HMMS_OUTPUT_TYPE] [--validation]
                    [-p] [-db DATABASE] [-s SNAKEFILE] [-t THREADS]
                    [-hm HMM_MODELS] [--concat_hmm_models] [--unlock]
                    [-w WORKFLOW] [-c CONFIG_FILE] [-v]

M-PARTY's main script

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        input FASTA file containing a list of protein
                        sequences to be analysed
  --input_seqs_db_const INPUT_SEQS_DB_CONST
                        input a FASTA file with a set of sequences from which
                        the user wants to create the HMM database from scratch
  -db DATABASE, --database DATABASE
                        FASTA database to run against the also user inputted
                        sequences. DIAMOND is performed in order to expand the
                        data and build the models. PlastEDMA has no in-built
                        database for this matter. If flag is given, download
                        of the default database will start and model built
                        from that. Defaults to UniProt DataBase.
  --hmm_db_name HMM_DB_NAME
                        name to be assigned to the hmm database to be created.
                        Its recomended to give a name that that describes the
                        family or other characteristic of the given sequences
  -it INPUT_TYPE, --input_type INPUT_TYPE
                        specifies the nature of the sequences in the input
                        file between 'protein', 'nucleic' or 'metagenome'.
                        Defaults to 'protein'
  -o OUTPUT, --output OUTPUT
                        name for the output directory. Defaults to
                        'PlastEDMA_results'
  --output_type OUTPUT_TYPE
                        choose report table outpt format from 'tsv', 'csv' or
                        'excel'. Defaults to 'tsv'
  -rt, --report_text    decides whether to produce or not a friendly report in
                        txt format with easy to read information
  --hmms_output_type HMMS_OUTPUT_TYPE
                        chose output type of hmmsearch run from 'out', 'tsv'
                        or 'pfam' format. Defaults to 'tsv'
  --validation          decides whether to perform models validation and
                        filtration with the 'leave-one-out' cross validation
                        methods. Call to set to True. Defaults to False
  -p, --produce_inter_tables
                        call if user wants to save intermediate tables as
                        parseale .csv files (tables from hmmsearch results
                        processing)
  --negative_db NEGATIVE_DB
                        path to a user defined negative control database.
                        Default use of human gut microbiome
  -s SNAKEFILE, --snakefile SNAKEFILE
                        user defined snakemake workflow Snakefile. Defaults to
                        '/workflow/Snakefile
  -t THREADS, --threads THREADS
                        number of threads for Snakemake to use. Defaults to 1
  -hm HMM_MODELS, --hmm_models HMM_MODELS
                        path to a directory containing HMM models previously
                        created by the user. By default PlastEDMA uses the
                        built-in HMMs from database in
                        'resources/Data/HMMs/After_tcoffee_UPI/'
  --concat_hmm_models   concatenate HMM models into a single file
  --unlock              could be required after forced workflow termination
  -w WORKFLOW, --workflow WORKFLOW
                        defines the workflow to follow, between "annotation",
                        "database_construction" and "both". Latter keyword
                        makes the database construction first and posterior
                        annotation. Defaults to "annotation"
  -c CONFIG_FILE, --config_file CONFIG_FILE
                        user defined config file. Only recommended for
                        advanced users. Defaults to 'config.yaml'. If given,
                        overrides config file construction from input
  --display_config      declare to output the written config file together
                        with results. Useful in case of debug
  -v, --version         show program's version number and exit
```