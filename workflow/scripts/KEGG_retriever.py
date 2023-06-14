from bs4 import BeautifulSoup
import requests
from requests.exceptions import HTTPError
import re
from tqdm import tqdm


def find_between(string, first, last):
    try:
        start = string.index(first) + len(first)
        end = string.index(last, start)
        return string[start:end]
    except ValueError:
        return ""

def get_kegg_genes(filepath: str, type_seq: str = "AA", ec_number = None, ko = None, verbose: bool = False):
    if ko == None and ec_number == None:
        raise ValueError("Either an E.C. number or a KO must be given")
    elif ko != None and ec_number != None:
        raise ValueError("Only give eith a E.C. number or KO")
    elif ko:
        ko_urlist = []
        for k in ko:
            ko_urlist.append(f'https://rest.kegg.jp/get/{k}')

            if verbose:
                print(f'Downloading {k} orthology genes')
        for link in ko_urlist:
            get_kegg_kosequences(filepath, link, korec = "ko", type_seq = type_seq, verbose=verbose)
    elif ec_number:
        ec_urlist = []
        for ec in ec_number:
            ec_urlist.append(f'https://rest.kegg.jp/get/{ec}')
            if verbose:
                print(f'Downloading {ec} gene')
        for link in ec_urlist:
            get_kegg_kosequences(filepath, link, korec = "ec", type_seq = type_seq, verbose=verbose)
    return filepath

def get_kegg_kosequences(filepath: str, url: str, korec: str = None, type_seq: str = "AA", verbose: bool = False):
    try:
        response = requests.get(url)
        soup = BeautifulSoup(response.text, "html.parser")
        genes = []
        start_line = False
        end_line = False
        with open(filepath, "a") as wf:
            for line in soup.get_text().split('\n'):
                if line.startswith("GENES"):
                    start_line = True
                    end_line = False
                    gene_line = line.split(None, 1)[1]
                    genes.extend(gene_line.split(';'))
                elif start_line and not end_line:
                    if line.startswith(' '):
                        gene_line = line.strip()
                        genes.extend(gene_line.split(';'))
                    else:
                        end_line = True
            # case of EC numbers with no acess to api genes
            if genes == []:
                if verbose:
                    print(f'Sequences from {url} not found. Trying in RefGene')
                wf.close()
                try:
                    get_kegg_refgene_sequences(filepath, url, korec, type_seq = type_seq, verbose = verbose)
                except:
                    print(f"[WARNING] Sequence for {url} not found")
            else:
                genes2 = []
                for entry in genes:
                    indiv = entry.split(" ")
                    indiv[0] = indiv[0].lower()
                    if len(indiv) == 2:
                        genes2.append("".join(indiv))
                    else:
                        for x in range(len(indiv)) + 1:
                            genes2.append(indiv[0] + indiv[x])
                print(genes2)
                for i in tqdm(range(len(genes)), desc = "Downloading genes"):
                    genes[i] = genes[i].replace(" ", "")
                    genes[i] = re.findall(".*:", genes[i])[0].lower() + re.findall(":.*", genes[i])[0][1:]
                    if type_seq == "AA":
                        url2 = f'https://rest.kegg.jp/get/{genes[i]}/aaseq'
                    else:
                        url2 = f'https://rest.kegg.jp/get/{genes[i]}/ntseq'
                    # if verbose:
                    #     print(f'Downloading sequence from {url2}')
                    response = requests.get(url2)
                    wf.write(response.text)
                if verbose:
                    print(f'Checking the existance of other genes for {url.split("/")[-1]} in RefSeq')
                get_kegg_refgene_sequences(filepath, url, korec, type_seq, verbose)
            wf.close()
    except:
        print(f"[WARNING] Sequence for {url} not found")
    
def get_kegg_refgene_sequences(filepath: str, url: str, korec: str = None, type_seq: str = "AA", verbose: bool = False):
    if korec == "ec":
        ec = url.split("/")[-1]
        url = f"https://www.genome.jp/dbget-bin/get_linkdb?-t+refgene+ec:{ec}"
    else:
        ko = url.split("/")[-1]
        print(ko)
        url = f"https://www.genome.jp/dbget-bin/get_linkdb?-t+refgene+ko:{ko}"
    try:
        response = requests.get(url)
    except HTTPError as http_err:
        print(f'HTTP error occurred: {http_err}')
    except Exception as err:
        print(f'Other error occurred: {err}')
    soup = BeautifulSoup(response.text, "html.parser")
    # Extract the gene IDs from the HTML
    gene_ids = [a.text for a in soup.find_all("a")]
    # Use the KEGG API to get the sequence information for a given list of gene IDs
    fasta = ""
    not_found = 0
    for gene_id in tqdm(gene_ids, desc = "Downloading genes", position=0, leave=True):
        if gene_id.startswith("RG"):
            # if verbose:
            #     print(f'Downloading {gene_id} gene')
            if type_seq == "AA":
                url = f"https://www.genome.jp/entry/-f+-n+a+{gene_id}"
            else:
                url = f"https://www.genome.jp/entry/-f+-n+n+{gene_id}"
            response = requests.get(url)
            text = response.text
            # Add the sequence information to the fasta string
            try:
                start = text.index("-->&gt;") + len("-->&gt;")
                end = text.index("</pre>", start)
                fasta += ">" + text[start:end]
            except:
                print(f"[WARNING] Sequence for {gene_id} not found")
                not_found += 1
                continue
        # print(fasta)
    file = open(filepath, "a")
    file.write(fasta)
    file.close()


# ec_number = ["3.1.1.101", "1.14.12.15", "1.18.1.-", "1.3.1.53"]
# ko = ["K18251", "K18252", "K18253", "K18254"]

# get_KEGG_sequences("/mnt/c/Users/Ze/Desktop/M-PARTY/.tests/KEGG_nuc_test.fasta", type="nuc", ko=ko, verbose=True)
# get_kegg_genes("/mnt/c/Users/Ze/Desktop/M-PARTY/.tests/KEGG_new_test.fasta", ec_number=ec_number, verbose = True)