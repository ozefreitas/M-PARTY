import requests
from requests.exceptions import HTTPError
import re
import time
from tqdm import tqdm
from mparty_util import get_soup


def get_IP_sequences(filepath: str, interpro_ID: str = None, reviewed: bool = False, protein: list = [], verbose: bool = False) -> str:
    """Function that will retrieve the protein sequences from InterPro if a compatible ID is given.

    Args:
        filepath (str): path to the file where the output will be writen
        interpro_ID (str, optional): ID for InterPro family. Defaults to None.
        reviewed (bool, optional): decides wheter to retrieve only reviewed sequences. Defaults to False.
        protein (list, optional): list of protein IDs of InterPro. Defaults to [].
        verbose (bool, optional): decides to print more information. Defaults to False.

    Raises:
        ValueError: first ValueError is raised in case neither parameter is given
        ValueError: second ValueError is raised in case neither parameter is given

    Returns:
        str: the path for the output file
    """
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
                if response.status_code == 408:
                    time.sleep(61)
                response.raise_for_status()
                json_resp = response.json()
            except HTTPError as http_err:
                print(f'HTTP error occurred: {http_err}')
            except Exception as err:
                print(f'Other error occurred: {err}')
            try:
                name = f'>{p}|{json_resp["metadata"]["name"]}|{json_resp["metadata"]["source_organism"]["scientificName"]}'
            except:
                name = p
            seq = json_resp["metadata"]["sequence"]
            list_prot[name] = seq
    elif interpro_ID != None:
        try:
            url = f'https://www.ebi.ac.uk/interpro/api/protein/UniProt/entry/InterPro/{interpro_ID[0]}/?extra_fields=sequence&page_size=100'
            response = requests.get(url)
            if response.status_code == 408:
                time.sleep(61)
            response.raise_for_status()
            json_resp = response.json()
        except HTTPError as http_err:
            print(f'HTTP error occurred: {http_err}')
        except Exception as err:
            print(f'Other error occurred: {err}')
        number = json_resp["count"]
        if verbose:
            print(f'Found InterPro entry for {interpro_ID} with {number} protein sequences')
            time.sleep(1)
        list_prot = {}
        counter = 0
        pbar = tqdm(total=number)
        while url:
            response = requests.get(url)
            if response.status_code == 408:
                time.sleep(61)
                continue
            response.raise_for_status()
            json_resp = response.json()
            for data in json_resp["results"]:
                if reviewed:
                    if data["metadata"]["source_database"] == "unreviewed":
                        continue
                seq = ""
                name4fasta = f'>{data["metadata"]["accession"]}|{data["metadata"]["name"]}|{data["metadata"]["source_organism"]["scientificName"]}'
                if verbose:
                    name = data["metadata"]["accession"]
                    print(f'Downloading {name}')
                for entry in data["entries"]:
                    for locations in entry["entry_protein_locations"]:
                        for frags in locations["fragments"]:
                            ini = frags["start"]
                            fin = frags["end"]
                            seq += data["extra_fields"]["sequence"][ini:fin]
                list_prot[name4fasta] = seq
                counter += 1
            url = json_resp["next"]
            pbar.update(1)
        pbar.close()
    with open(filepath, "w") as wf:
        for k, v in list_prot.items():
            wf.write(k + "\n")
            fasta_seq = ""
            for i in range(0, len(v), 60):
                # print(i)
                if i + 60 > len(v):
                    fasta_seq += v[i:] + "\n"
                    # print(fasta_seq)
                    break
                fasta_seq += v[i:i+60] + "\n"
            # print(fasta_seq)
            wf.write(fasta_seq)
    wf.close()
    return filepath


# get_gene_IDS(interpro_ID = "ola", protein = ["A0A000", "A0A001", "A0A002"])
# print(get_gene_IDS(interpro_ID = "IPR000006", verbose=True))
# print(get_gene_IDS(interpro_ID = "IPR000003"))
# dicionario = get_gene_IDS(protein = ["A0A000", "A0A001", "A0A002"])
# create_fasta(dicionario, "/mnt/c/Users/Ze/Desktop/M-PARTY/.tests/test_interpro.fasta")
# get_IP_sequences("/mnt/c/Users/Ze/Desktop/M-PARTY/.tests/test_interpro.fasta", protein = ["A0A000", "A0A001", "A0A002"], verbose=True)