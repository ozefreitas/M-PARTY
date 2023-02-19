from bs4 import BeautifulSoup
import requests
import re


def find_between(string, first, last):
    try:
        start = string.index(first) + len(first)
        end = string.index(last, start)
        return string[start:end]
    except ValueError:
        return ""


def get_gene_ids(ec_number: str = None, ko: str = None) -> list:
    # Use the KEGG API to get the gene IDs for a given KEGG entry
    if ko == None and ec_number == None:
        raise ValueError("Either an E.C. number or a KO must be given")
    elif ko != None and ec_number != None:
        raise ValueError("Only give eith a E.C. number or KO")
    elif ec_number:
        url = f"https://www.genome.jp/dbget-bin/get_linkdb?-t+refgene+ec:{ec_number}"
    elif ko:
        url = f"https://www.genome.jp/dbget-bin/get_linkdb?-t+refgene+ko:{ko}"
    response = requests.get(url)
    soup = BeautifulSoup(response.text, "html.parser")
    # Extract the gene IDs from the HTML
    gene_ids = [a.text for a in soup.find_all("a")]
    return gene_ids


def get_gene_sequences(gene_ids):
    # Use the KEGG API to get the sequence information for a given list of gene IDs
    fasta = ""
    not_found = 0
    for gene_id in gene_ids:
        if gene_id.startswith("RG"):
            # print(gene_id)
            # url = f'https://www.rest.kegg.jp/get/{gene_id}/aaseq'
            url = f"https://www.genome.jp/entry/-f+-n+a+{gene_id}"
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
    return fasta, not_found


def save_fasta(filepath: str, fasta: str):
    file = open(filepath, "w")
    file.write(fasta)
    file.close()


ec_number = "3.1.1.101"
ko = "K21104"
ref_genes = get_gene_ids(ec_number = ec_number)
fasta, not_found = get_gene_sequences(ref_genes)
save_fasta("/mnt/c/Users/Ze/Desktop/M-PARTY/.tests/KEGG_test.fasta", fasta)

ref_genes_KO = get_gene_ids(ko=ko)
fasta_KO, not_found = get_gene_sequences(ref_genes_KO)
save_fasta("/mnt/c/Users/Ze/Desktop/M-PARTY/.tests/KEGG_test_KO.fasta", fasta_KO)