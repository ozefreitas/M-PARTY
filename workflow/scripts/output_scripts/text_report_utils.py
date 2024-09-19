def write_var_file():
    variables = {"number_init_hmms": 0, 
            "query_names": "", 
            "number_validated_hmms": 0, 
            "bit_threshold": 0, 
            "eval_threshold": 0, 
            "inputed_seqs": 0, 
            "unique_seqs": 0}
    return variables


def write_text_report(config, path, args, variables):

    if config["hmm_validation"] == True:
        with open(path + "text_report.txt", "w") as f:
            f.write(f"M-PARTY hits report:\n \
                    \nFrom a total number of {number_init_hmms} HMM profiles initially considered, only {len(query_names)} where considered"
                    f"for the final report.\nUser defined validation to true with {args.negative_db} database, from which resulted"
                    f"in {number_validated_hmms}. After annotation, another filtering process was performed considering the values from bit score and E-value from the HMM search run, \n"
                    f"in which the considered bit score threshold was {bit_threshold} and E-value was {eval_threshold}.\n"
                    f"Also, {len(inputed_seqs)} initial query sequences where inputed form {args.input} file, from which {len(unique_seqs)}\n"
                    f"out of these {len(inputed_seqs)} were considered to have a hit against the HMM database.")
            f.close()
    else:
        with open(path + "text_report.txt", "w") as f:
            f.write(f"M-PARTY hits report:\n \
                    \nFrom a total number of {number_init_hmms} HMM profiles initially considered, only {len(query_names)} where considered"
                    "for the final report.\nFiltering process was performed considering the values from bit score and E-value from the HMM search run, \n"
                    f"in which the considered bit score threshold was {bit_threshold} and E-value was {eval_threshold}.\n"
                    f"Also, {len(inputed_seqs)} initial query sequences where inputed form {args.input} file, from which {len(unique_seqs)}\n"
                    f"out of these {len(inputed_seqs)} were considered to have a hit against the HMM database.")
            f.close()