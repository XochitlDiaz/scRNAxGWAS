# 16 Feb 2024
# Transform gene sympols to HGNC if necessary
# Xochitl Diaz


# Fix later
import requests
from bioservices import BioMart
import os
import csv

def check_hgnc_gene(gene_name):
    url = f'https://www.genenames.org/cgi-bin/download/custom?col=gd_app_sym&col=gd_hgnc_id&status=Approved&hgnc_dbtag=on&order_by=gd_app_sym_sort&format=text&submit=submit&custom_download=1&filename=&gz=0&limit=&hgnc_dbtag=on&multiple=1&query={gene_name}'
    response = requests.get(url)
    return response.ok and len(response.text.strip().split('\n')) > 1
#https://www.genenames.org/data/gene-symbol-report/#!/hgnc_id/HGNC:6149


### TEST ###
genes = pd.read_csv("./dummy_data/mic_genes.txt", header= None)
genes=list(genes[0])
for gene in genes:
    print(gene)
# https://www.genenames.org/tools/multi-symbol-checker/
    

def convert_gene_names_to_ensembl_ids(gene_names):
    biomart = BioMart(verbose=False)
    biomart.datasets(verbose=False)
    
    # Select the appropriate dataset (e.g., human genes)
    biomart.dataset = 'hsapiens_gene_ensembl'
    
    # Define the attributes you want to retrieve (ENSEMBL Gene ID and HGNC Symbol)
    attributes = ['ensembl_gene_id', 'external_gene_name']
    
    # Query BioMart to retrieve ENSEMBL gene IDs corresponding to the given gene names
    results = biomart.query(attributes, filters={'external_gene_name': gene_names})
    
    # Organize the results into a dictionary where gene names map to ENSEMBL gene IDs
    gene_names_to_ensembl = {}
    for row in results:
        ensembl_id, gene_name = row
        gene_names_to_ensembl[gene_name] = ensembl_id
    
    return gene_names_to_ensembl

def convert_ensembl_to_hgnc(ensembl_ids):
    biomart = BioMart(verbose=False)
    biomart.datasets(verbose=False)
    
    # Select the appropriate dataset (e.g., human genes)
    biomart.dataset = 'hsapiens_gene_ensembl'
    
    # Define the attributes you want to retrieve (ENSEMBL Gene ID and HGNC Symbol)
    attributes = ['ensembl_gene_id', 'hgnc_symbol']
    
    # Query BioMart to retrieve HGNC symbols corresponding to the given ENSEMBL IDs
    results = biomart.query(attributes, filters={'ensembl_gene_id': ensembl_ids})
    
    # Organize the results into a dictionary where ENSEMBL IDs map to HGNC symbols
    ensembl_to_hgnc = {}
    for row in results:
        ensembl_id, hgnc_symbol = row
        ensembl_to_hgnc[ensembl_id] = hgnc_symbol
    
    return ensembl_to_hgnc


#### TEST ####
# Example gene names
gene_file = open("./mic_genes.csv", "r")
gene_names = list(csv.reader(gene_file, delimiter=","))
gene_file.close()
# Convert gene names to ENSEMBL IDs
gene_names_to_ensembl = convert_gene_names_to_ensembl_ids(gene_names)

# Get the ENSEMBL IDs
ensembl_ids = list(gene_names_to_ensembl.values())

# Convert ENSEMBL IDs to HGNC symbols
ensembl_to_hgnc = convert_ensembl_to_hgnc(ensembl_ids)

# Print the results
for gene_name, ensembl_id in gene_names_to_ensembl.items():
    hgnc_symbol = ensembl_to_hgnc.get(ensembl_id, "Not found")
    print(f"{gene_name} ({ensembl_id}) -> {hgnc_symbol}")
