# Generate and standarize cell types across datasets to be compared
# 31 Jan 2024
# XochitlDiaz

# python 3.11.5



# NOTE: this celltype integration was previously made using Seurat software
# inside of R. The code to do this can be found on the cluster as: 

# in future verisons, said process should be included in this preprocessing code.

### GENRATE A FUNCTION TO BE IMPORTED INT OTHER PROCESSES TO 
### STANDARIZE THE DATA SETS THAT WILL BE USED

### DATA QUALITY CONTROL TO REMOVE LOW QUALITY CELLS

### IN PUT : RAW COUNTS MATRIX AND META DATA MATRIX, ORDERED GENE CODE ARRAY
### OUTPUT : WELL FORMATED ANDATA OBJECTS SAVED IN A FILE 
###          + MARKDOWN FILE OF QC MADE?


import scanpy as sc
from collections import Counter
import numpy as np
import pandas as pd
import scipy.sparse  
import anndata
import os
# If converting a .rds to anndata (not recommended)
import anndata2ri

# Addons from the UCLA tutorial
import matplotlib.pyplot as plt
import seaborn as sns

#If data is in csv
exp_mat = pd.read_csv("./dummy_data/counts.csv", header=None)

exp_mat= anndata.AnnData(exp_mat)

#If data is in .mtx
#ADsyn = scipy.io.mmread("../../../ADsyn/filtered_count_matrix.mtx")
#ADsyn = ADsyn.T
#sparse_coo = scipy.sparse.coo_matrix(ADsyn)  
#ADsyn = sparse_coo.tocsr()
#ADsyn = anndata.AnnData(ADsyn)

# Open var names file and assign to anndata
gene_names= open("./dummy_data/gene_names.txt").read().split("\n")
exp_mat.var_names = gene_names


# Open metadata file and assign to anndata
metad= pd.read_csv("./dummy_data/metadata.csv")
exp_mat.obs= metad


# Write .h5ad object
f5ad_file_name = input("Data name:")
f5ad_file_name = f5ad_file_name + ".h5ad"
exp_mat.write(f5ad_file_name)