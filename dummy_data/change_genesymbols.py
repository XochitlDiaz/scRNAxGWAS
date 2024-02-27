# 16 Feb 2024
# transform gene sympols to HGNC if necessary
# Xochitl Diaz


# Fix later
import requests


def check_hgnc_gene(gene_name):
    url = f'https://www.genenames.org/cgi-bin/download/custom?col=gd_app_sym&col=gd_hgnc_id&status=Approved&hgnc_dbtag=on&order_by=gd_app_sym_sort&format=text&submit=submit&custom_download=1&filename=&gz=0&limit=&hgnc_dbtag=on&multiple=1&query={gene_name}'
    response = requests.get(url)
    return response.ok and len(response.text.strip().split('\n')) > 1
https://www.genenames.org/data/gene-symbol-report/#!/hgnc_id/HGNC:6149


### TEST ###
genes = pd.read_csv("./dummy_data/gene_names.txt", header= None)
genes=list(genes[0])
for gene in genes:
    print(gene)
# https://www.genenames.org/tools/multi-symbol-checker/
    
