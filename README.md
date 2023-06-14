# M-PARTY - Mining Protein dAtabases foR Targueted enzYmes
#### Version 1.0.0


<br>
M-PARTY is a free to use, open source user friendly CLI (early release) implemented workflow and database for the detection of plastic degrading enzymes in metagenomic samples, through structural annotation using Hidden Markov Models.
<p>
<br>

## Index

1. [Introduction](https://github.com/ozefreitas/M-PARTY#introduction)
2. [Installation](https://github.com/ozefreitas/M-PARTY#installation)
    1. [GitHub](https://github.com/ozefreitas/M-PARTY#github)
    2. [Bioconda](https://github.com/ozefreitas/M-PARTY#bioconda)
3. [Usage](https://github.com/ozefreitas/M-PARTY#usage)
    1. [Annotation](https://github.com/ozefreitas/M-PARTY#annotation-workflow)
    2. [Database Construction](https://github.com/ozefreitas/M-PARTY#database-construction)
    3. [Validation](https://github.com/ozefreitas/M-PARTY#validation)
    4. [Full M-PARTY Execution](https://github.com/ozefreitas/M-PARTY#full-mparty-execution)
    5. [Metagenomic Analisys](https://github.com/ozefreitas/M-PARTY#metagenomic-analisys)
4. [Output](https://github.com/ozefreitas/M-PARTY#output)
5. [Additional arguments](https://github.com/ozefreitas/M-PARTY#additional-arguments)

<br>

## Introduction 

M-PARTY is a free to use, open source user friendly CLI implemented workflow and database for the detection of plastic degrading enzymes in metagenomic samples, through structural annotation using Hidden Markov Models, that allows the user to freely interacte with the tool in-built databases and backbone. <p>
M-PARTY compreends a extensive HMM database, built with state of the art checked enzymatic sequences able to degrade plastic polymers, which is used to carry the structural annotation of given sequences. <p>
First version of M-PARTY is only available for mining PE (polyethylene), as latter version will compreend a more vast list of plastics to analyse. <p>
Also, M-PARTY is meant to analyse metagenomic sequences, but version 0.1.0 only accepts single FASTA aminoacidic sequences. Basic steps of M-PARTY annotation workflow in its frist stages are: 

1. The acceptence of any number of protein sequences in a single FASTA file as query; also an KEGG ID representing the protein 
sequences involved in certain reaction; an InterPro ID or a set of protein IDs of interest from latter database;
2. Execution of hmmsearch from the [HMMER package](https://www.hmmer.org/) using the built HMMs from previously knowns sequences able to have some kind of PE deterioration levels as database; and [KMA](https://bitbucket.org/genomicepidemiology/kma/src/master/) to map and search raw metagenomes for interest genes; 
3. Able to input three different genres of files:
    - Protein datasets
    - Assembled metagenomes
    - Raw metagenomes
4. A quality benchmark to determine good and bad hits from the queries against models;
5. Three output files, consisting in a FASTA file with the protein sequences returned as a hit from the search, a report in text format (if requested by the user) with simply puted information about the inputed and already built (HMMs) data, run and processing parameters and conclusions, and an easy to read report table in xlsx format, with all the important data about the annotation results, in particular:
    - Sequence IDs
    - HMM IDs (Degraded plastic + number)
    - Bit scores
    - E-values

<br>

## Installation
<p>

### GitHuh Clonning 
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

<br>

### Bioconda
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
<p>

#### Annotation workflow

The main and most basic use for M-PARTY is the annotation with Hidden Markov Models (must be previously created by the user):<p>

```
m-party.py -i path/to/input_file -o output_folder -rt --output_type excel --hmm_db_name <db_name> --verbose
```

where the **`-i` input file** must be in FASTA format and contain only (for the time being) aminoacidic sequences, otherwise, program will exit. **`-o` output folder** can be a pre-existing folder or any name for a folder that will be created anyways. The **`-rt`** option flag instructs the tool to include in the output the report in text format, for an easier interpretation of the annotation results and conclusion taking. Also, **`--output_type`** is recommended to be set to "excel" on these earlier versions, as other output format for the table report will be incrementally coded. **`--hmm_db_name`** is mandatory and represents a name to be given to each run, and folders with that name will be saved as databases.

<br>

### Database Construction <p>

M-PARTY does not have any pre-built database, and so, all HMMs must be generated from scratch from a given set of proteins/nucleotides. M-PARTY accepts this input from 3 distinct methods:

  1. A FASTA file with sequences with known functions from the user;
  2. A KEGG Orthodology ID(s) (KO) or E.C. number(s);
  3. An InterPro ID or Protein ID(s)

With previous knowledge of a given reaction or protein family, a user can input an ID which represents a set of sequences involved in such reaction, from where this tool will automatically search and download.

<br>

**FASTA file**

If you possess a FASTA file with interest sequences for your study to after be searched:

```
m-party.py -w database_construction --input_seqs_db_const path/to/interest_sequences --hmm_db_name <db_name>
```
<br>

**KEGG**

If you want to build an HMM database from a reaction represented by a certain KO:

```
m-party.py -w database_construction --kegg <KO> --hmm_db_name <db_name>
``` 
or for a certain E.C. number (EC):

```
m-party.py -w database_construction --kegg <EC> --hmm_db_name <db_name>
``` 

By default, M-PARTY will download the aminoacid sequences of the found enzymes for each ID. If you wish to build the models for nucleic sequences just add the **`--input_type_db_const`** and set it to "nucleic" like (in both cases of FASTA user file or KEGG IDs):

```
m-party.py -w database_construction --input_seqs_db_const path/to/interest_sequences OR --kegg <EC/KO> --input_type_db_const nucleic --hmm_db_name <db_name> 
``` 
<br>

**InterPro**

As for the case of InterPro retriver, just change the argument kegg for **`--interpro´**:

```
m-party.py -w database_construction --interpro <IPR> --hmm_db_name <db_name>
```
or
```
m-party.py -w database_construction --interpro <PID> --hmm_db_name <db_name>
```

For the case of an InterPro ID (IPR) or a Protein ID (PID). In this database, only aminoacid sequences are available, and so it is not possible to set **`--input_type_db_const`** to "nucleic".

Interpro its not a curated database and so most of his entries are unreviewed. To counter this, you can add the **`--curated`** flag to the previous line

*Note:* All commands must have the **`--hmm_db_name <db_name>`** argument! Otherwise, M-PARTY will instantly raise a `ValueERROR`.

<br>

### Validation <p>
M-PARTY also has available a validation workflow, as some models can be false positive and give deceiving results. So giving the name where the HMM are stored:
```
m-party.py --hmm_validation --hmm_db_name <db_name>
```

 Obviously, using the 'leave-one-out' cross validation method, a dataset of negative sequences must be given, and so each dataset will be distinct in each run, depending of the content of the HMMs. So just add the **`--negative_db`** with a FASTA dataset of sequences that you know to be different from the ones of interest:

 ```
m-party.py --hmm_validation --hmm_db_name <db_name> --negative_db path/to/negative_dataset
```

Now, if you want, you can instantly run the annotation workflow from a set of proteins of your liking, and so performing the validation beforehand, only if you already ran the *database construction* workflow, and so the models are already present, with:

```
m-party -i path/to/input_file -o path/to/output_folder -rt --output_type excel --hmm_db_name <db_name> --hmm_validation 
```

**Full MPARTY Execution** <p>

How could not miss, all this can be done with a single command, from the construction of the models, to the annotation of the unknown sequence file.

```
m-party.py -w both -i path/to/input_file (--input_seqs_db_const path/to/interest_sequences OR --kegg <EC/KO> OR --interpro <IPR>) -o path/to/output_folder -rt --output_type excel --hmm_db_name <db_name> --verbose 
```

### Metagenomic Analisys <p>

M-PARTY also suports the search of genes in metagenome samples:

```
m-party.py -w database_construction --it metagenome --kegg <KO> --input_type_db_const nucleic --hmm_db_name <db_name>
```

I firstly recommend to run the *database construction* workflow in order to download or input the interest genes. In order to avoid building the models and so waste time, **`--it`** is set to metagenome.<p>
After this, just give the metagenome file, say again that what is inside is in fact a metagenome, give the same <db_name> and run it:

```
m-party.py -i path/to/metagenome -o path/to/output_folder -it metagenome --hmm_db_name <db_name>
```
*Warning:* This method is only viable for nucleotide sequences, so **`--input_type_db_const nucleic`** is obligatory!

Other way is to just do both workflows at the same time:

```
m-party.py -w both -it metagenome -i path/to/metagenome -o path/to/output_folder --kegg <KO> --input_type_db_const nucleic --hmm_db_name <db_name> --verbose
```

Note that instead of **`--kegg`** you can put instead the **`--input_seqs_db_const`** argument with a FASTA file of nucleotides.

<br>

## Output
<p>

M-PARTY will result in three distinct outputs: **report table**, **text report** and **aligned**. In earlier versions, **report table** is only available in *excel* format, although later will also be for *tsv* and *csv*. **Text report** is a user friendly easy to understand summary of the annotation run performed by M-PARTY, and embrace a series of useful information for the user, depending on the given arguments. For last, **aligned** is a FASTA file with all the sequences that had a match in one or more models (this will be refined as model benchmarking and validation are introduced into M-PARTY). 

Optional output file is the *config file*, a file with all parameters used in each M-PARTY run, and that can be useful to traceback error. If you are having trouble generating the expect content in the other output files, you can add the **`--diplay_config`** flag to every command you execute.

<br>

## Aditional arguments
<p>

M-PARTY is in continuous development, so the **"validation"** workflow is needing some changes, as well as the expansion argument, which still needs to be reviewed and validated. So I highly recommend you to follow the steps in the [usage](https://github.com/ozefreitas/M-PARTY#usage) section.

```
usage: m-party.py [-h] [-i INPUT] [--input_seqs_db_const INPUT_SEQS_DB_CONST]
                  [-db DATABASE] [--hmm_db_name HMM_DB_NAME] [-it INPUT_TYPE]
                  [--input_type_db_const INPUT_TYPE_DB_CONST] [--consensus]
                  [-o OUTPUT] [--output_type OUTPUT_TYPE] [-rt]
                  [--hmms_output_type HMMS_OUTPUT_TYPE] [--hmm_validation]
                  [-p] [--negative_db NEGATIVE_DB] [-s SNAKEFILE] [-ex]
                  [--kegg KEGG [KEGG ...]]
                  [--interpro INTERPRO [INTERPRO ...]] [--curated]
                  [-t THREADS] [--align_method ALIGN_METHOD]
                  [--aligner ALIGNER] [-hm HMM_MODELS] [--concat_hmm_models]
                  [--unlock] [-w WORKFLOW] [-c CONFIG_FILE] [--overwrite]
                  [--verbose] [--display_config] [-v]

M-PARTY's main script

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        input FASTA file containing a list of protein
                        sequences to be analysed
  --input_seqs_db_const INPUT_SEQS_DB_CONST
                        input a FASTA file with a set of sequences from which
                        the user wants to create the HMM database from
                        scratch.
  -db DATABASE, --database DATABASE
                        FASTA database to run against the also user inputted
                        sequences. DIAMOND is performed in order to expand the
                        data and build the models. M-PARTY has no in-built
                        database for this matter. If flag is given, download
                        of the default database will start and model built
                        from that. Defaults to UniProt DataBase.
  --hmm_db_name HMM_DB_NAME
                        name to be assigned to the hmm database to be created.
                        Its recomended to give a name that that describes the
                        family or other characteristic of the given sequences.
                        Be carefull as what name to use, as this will define
                        what HMMs will be used for the search
  -it INPUT_TYPE, --input_type INPUT_TYPE
                        specifies the nature of the sequences in the input
                        file between 'protein', 'nucleic' or 'metagenome'.
                        Defaults to 'protein'
  --input_type_db_const INPUT_TYPE_DB_CONST
                        specifies the nature of the input sequences for the
                        database construction between 'nucleic' and 'protein'.
                        Defaults to 'protein'.
  --consensus           call to build consensus sequences when building the
                        database, in order to run KMA against raw metagenomes
  -o OUTPUT, --output OUTPUT
                        name for the output directory. Defaults to
                        'MPARTY_results'
  --output_type OUTPUT_TYPE
                        choose report table outpt format from 'tsv', 'csv' or
                        'excel'. Defaults to 'tsv'
  -rt, --report_text    decides whether to produce or not a friendly report in
                        txt format with easy to read information
  --hmms_output_type HMMS_OUTPUT_TYPE
                        chose output type of hmmsearch run from 'out', 'tsv'
                        or 'pfam' format. Defaults to 'tsv'
  --hmm_validation      decides whether to perform models validation and
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
  -ex, --expansion      Decides wheter to expand the interest dataset.
                        Defaults to False.
  --kegg KEGG [KEGG ...]
                        input KEGG ID(s) to download respective sequences, in
                        order to build a pHMM based on those
  --interpro INTERPRO [INTERPRO ...]
                        input InterPro ID(s) to download the respective
                        sequences, in order to build a pHMM based on those
  --curated             call to only retrieve reviewed sequences from InterPro
  -t THREADS, --threads THREADS
                        number of threads for Snakemake to use. Defaults to
                        max number of available logical CPUs.
  --align_method ALIGN_METHOD
                        chose the alignment method for the initial sequences
                        database expansion, between 'diamond', 'blast' and
                        'upimapi'. Defaults to 'upimapi'
  --aligner ALIGNER     chose the aligner program to perform the multiple
                        sequence alignment for the models between 'tcoffee'
                        and 'muscle'. Defaults to 'tcoffee'.
  -hm HMM_MODELS, --hmm_models HMM_MODELS
                        path to a directory containing HMM models previously
                        created by the user. By default M-PARTY uses the
                        built-in HMMs from database in
                        'resources/Data/HMMs/After_tcoffee_UPI/'
  --concat_hmm_models   call to not concatenate HMM models into a single file.
                        Defaults to True
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
  --overwrite           Call to overwrite inputted files. Defaults to False
  --verbose             Call so M-PARTY display more messaging
  --display_config      declare to output the written config file together
                        with results. Useful in case of debug
  -v, --version         show program's version number and exit
```