# Generate Cell type Programs
# 1 Feb 2024
# XochitlDiaz

# This function will generate Cell Type programs out 
# of a previously QCed .h5ad file

import scanpy as sc
from collections import Counter
import numpy as np
import pandas as pd
import scipy.sparse  
import anndata
import os


# In future versions make this into a a propper input
# .h5ad data file with celltype id and gene names
exp_mat_dir = "./dummy_data/simdata.h5ad"
# directory format to save the output files
save_dir = "./dummy_data/celltype_program"
# on the metadata part of .h5ad what is the column name for the annotated celltype
# Can be cheked by print(exp.obs.columns)
ct_colname = 'Celltype'
# Indicate the name of the column that defines diseare status
status_colname = 'Status'
scdatasets = 'Status'
# Indicate the name of the cloumn that defines patient ID
sample_colname = 'Patient'

def genCellTypePrograms(exp_mat_dir, ct_colname, status_colname, sample_colname):
    #read in data
    exp = sc.read_h5ad(exp_mat_dir)
    #change referenced column names for anndata manipulation
    exp.obs.rename(columns={ct_colname:'cell_type', status_colname:'status', sample_colname:'patient_id'},inplace=True)
    
    #######
    # When gnerating Celltype Data, first give some information
    # on the data that is being processed
    print("Generating celltype programs for the file:" + exp_mat_dir)
    data_name = exp_mat_dir[exp_mat_dir.rfind("/") + 1:(exp_mat_dir.rfind(".h5ad"))]
    print("Automatically chosen data name:  " + data_name)
    print("\nDescription of the data:")
    print(str(exp.shape[0]) + " cells | " + str(exp.shape[1]) + " genes | "+ str(len(set(exp.obs['patient_id']))) + " samples")
    print(str(len(set(exp.obs['cell_type']))) + " cell types: " + str(list(exp.obs['cell_type'].unique())))
    print("cell state groups: " + str(list(exp.obs['status'].unique())))
    ######

    delabel =  'celltype_DE'     
    # the data should have more than 10 cells of the cell type to study
    counts = Counter(exp.obs["cell_type"])
    #minimum cell type frequency
    min_ctfreq= 10
    frequent_ct = [key for key, count in counts.items() if count > 10]
    adata = exp[exp.obs['cell_type'].isin(frequent_ct)]

    #Differential expression analysis one celltype vs rest
    #Rank genes for characterizing groups
    sc.tl.rank_genes_groups(adata, groupby='cell_type', key_added='celltype_DE', use_raw=False, method='wilcoxon', n_genes= adata.shape[1])
    #Plot ranked genes
    sc.pl.rank_genes_groups(adata, n_genes=25, key='celltype_DE', sharey=False)
    exp.uns['celltype_DE'] = adata.uns['celltype_DE']
    
    #Write raw results of DE to file 
    genes = list(set(exp.var_names))
    gene2idx = {gene:i for i, gene in enumerate(genes)}
    pvalmtxs, logfoldmtxs, scoremtxs = [], [], []

    cellsubsets = exp.uns[delabel]['names'].dtype.fields.keys()

    # Which cells in the dataset are from which type
    cell2idx = {cellsubset: mask for mask, cellsubset in enumerate(cellsubsets)}

    # create empty matrix
    # adj pvalue of each gene for each group
    pvalmtx = np.zeros((len(gene2idx), len(cell2idx)))
    # logfold change of each gene for each group
    logfoldmtx = np.zeros((len(gene2idx), len(cell2idx)))
    # zscore of each gene for each group
    scoremtx = np.zeros((len(gene2idx), len(cell2idx)))

    # loop through and fill up the matrix with pvalue, logfold and score
    iters =  zip(exp.uns[delabel]['names'], 
                exp.uns[delabel]['pvals_adj'], 
                exp.uns[delabel]['logfoldchanges'], 
                exp.uns[delabel]['scores'])
    for gene, pval, logfold, score in iters:
        for cell_subset in cellsubsets:
            if gene[cell_subset] in gene2idx:
                pvalmtx[gene2idx[gene[cell_subset]], cell2idx[cell_subset]] = pval[cell_subset]
                logfoldmtx[gene2idx[gene[cell_subset]], cell2idx[cell_subset]] = logfold[cell_subset]
                scoremtx[gene2idx[gene[cell_subset]], cell2idx[cell_subset]] = score[cell_subset]

    # matrix to dataframe
    cellsubsets = [ct+"_celltype" for ct in cellsubsets]
    pvalmtxs.append(pd.DataFrame(pvalmtx, index=genes, columns=cellsubsets))
    logfoldmtxs.append(pd.DataFrame(logfoldmtx, index=genes, columns=cellsubsets))
    scoremtxs.append(pd.DataFrame(scoremtx, index=genes, columns=cellsubsets))
    
    pvalmtxs = pd.concat(pvalmtxs, axis=1)
    logfoldmtxs = pd.concat(logfoldmtxs, axis=1)
    scoremtxs = pd.concat(scoremtxs, axis=1)

    # write matrix to file
    # This matrix contains a cloumn per cell type and rows per genes
    # It stores the results of DE of one celltype vs all other celltypes
    os.mkdir("./celltype_program")
    save_dir="./celltype_program"
    pvalmtxs.to_csv("%s/%s_adjpval.csv"%(save_dir, data_name))
    logfoldmtxs.to_csv("%s/%s_logfold.csv"%(save_dir, data_name))
    scoremtxs.to_csv("%s/%s_zscore.csv"%(save_dir, data_name))

    # Transform p-value to x=-2log(P)
    scoremtxs = scoremtxs.clip(lower=0)
    scoremtxs = pd.DataFrame(scipy.stats.norm.sf(scoremtxs), index=scoremtxs.index, columns=scoremtxs.columns)
    scoremtxs = scoremtxs+1e-08
    scoremtxs = -2*np.log(scoremtxs)

    # Minimization
    scoremtxs = (scoremtxs - scoremtxs.min())/(scoremtxs.max()-scoremtxs.min())

    #Save transformed scores
    scoremtxs.to_csv("%s/%s_transgenescores.csv"%(save_dir, data_name))

