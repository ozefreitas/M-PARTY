import requests
from requests.exceptions import HTTPError
import re


def get_gene_IDS(interpro_ID: str = None, protein: list = [], verbose: bool = False):
    if interpro_ID == None and protein == []:
        raise ValueError("Either an InterPro ID or a list with protein IDs must be given")
    elif interpro_ID != None and protein != []:
        raise ValueError("Only give an InterPro ID or a list with protein IDs")
    elif protein != []:
        list_prot = {}
        for p in protein:
            try:
                url = f'https://www.ebi.ac.uk/interpro/api/protein/UniProt/{p}'
                response = requests.get(url)
                response.raise_for_status()
                json_resp = response.json()
            except HTTPError as http_err:
                print(f'HTTP error occurred: {http_err}')
            except Exception as err:
                print(f'Other error occurred: {err}')
            try:
                name = ">" + p + "|" + json_resp["metadata"]["name"] + "|" + json_resp["metadata"]["source_organism"]["scientificName"]
            except:
                name = p
            seq = json_resp["metadata"]["sequence"]
            list_prot[name] = seq
        return list_prot
    elif interpro_ID != None:
        try:
            url = f'https://www.ebi.ac.uk/interpro/api/entry/InterPro/{interpro_ID}'
            response = requests.get(url)
            response.raise_for_status()
            json_resp = response.json()
        except HTTPError as http_err:
            print(f'HTTP error occurred: {http_err}')
        except Exception as err:
            print(f'Other error occurred: {err}')
        if verbose:
            print(f'Found InterPro entry for {json_resp["metadata"]["name"]["name"]}')
        try:
            url = f'https://www.ebi.ac.uk/interpro/api/entry/interpro/{interpro_ID}/protein/uniprot'
            response = requests.get(url)
            response.raise_for_status()
            json_resp = response.json()           
        except HTTPError as http_err:
            print(f'HTTP error occurred: {http_err}')
        except Exception as err:
            print(f'Other error occurred: {err}')

        return len(json_resp["protein_subset"])


def create_fasta(seq_dict: dict, filename: str):
    with open(filename, "w") as wf:
        for k, v in seq_dict.items():
            wf.write(k)
            fasta_seq = ""
            for i in range(0, len(v), 60):
                if i + 60 > len(v):
                    fasta_seq += v[i:] + "\n"
                    print(fasta_seq)
                    break
                fasta_seq += v[i:i+60] + "\n"
            wf.write(fasta_seq)
    wf.close()


# get_gene_IDS(interpro_ID = "ola", protein = ["A0A000", "A0A001", "A0A002"])
# print(get_gene_IDS(interpro_ID = "IPR000003"))
dicionario = get_gene_IDS(protein = ["A0A000", "A0A001", "A0A002"])
create_fasta(dicionario, "/mnt/c/Users/Ze/Desktop/M-PARTY/.tests/test_interpro.fasta")