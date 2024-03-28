import scanpy as sc
from collections import Counter
import numpy as np
import pandas as pd
import scipy.sparse  
import anndata
import os


#If data is in csv
exp_mat = pd.read_csv("../testdata/mic_counts.csv", header=None)

exp_mat= anndata.AnnData(exp_mat)

# Open var names file and assign to anndata
gene_names= open("../testdata/mic_genes.csv").read().split("\n")
exp_mat.var_names = gene_names


# Open metadata file and assign to anndata
metad= pd.read_csv("../testdata/mic_metadata.csv")
exp_mat.obs= metad


# Write .h5ad object
f5ad_file_name = "mic_data"
f5ad_file_name = f5ad_file_name + ".h5ad"
exp_mat.write_h5ad("../testdata/" + f5ad_file_name)
