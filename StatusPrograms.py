# Generate Cell type Programs
# 7 Feb 2024
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
import math


# In future versions make this into a a propper input
# .h5ad data file with celltype id and gene names
exp_mat_dir = "./dummy_data/simdata.h5ad"
# directory format to save the output files
save_dir = "./dummy_data/status_program"
# on the metadata part of .h5ad what is the column name for the annotated celltype
# Can be cheked by print(exp.obs.columns)
ct_colname = 'Celltype'
# Indicate the name of the column that defines diseare status
status_colname = 'Status'
scdatasets = 'Status'
# Indicate the name of the cloumn that defines patient ID
sample_colname = 'Patient'
#In the Statues column, what is the string that idicates disease (opposite of reference for comparison)
diseaselab = 'Case'

def genStatusPrograms(exp_mat_dir, ct_colname, status_colname, sample_colname,diseaselab):
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
    deadatas = {}

    # Assuming there can be more than one disease state,
    # this loop will perform all possible paird comparisons,
    # according to cell type and a given disease state 
    states= list(exp.obs['status'].unique())
    base_state = next((s for s in states if s != diseaselab), None)

    # Set up matrix to write raw results of DE on
    genes = list(set(exp.var_names))
    gene2idx = {gene:i for i, gene in enumerate(genes)}
    pvalmtxs, logfoldmtxs, scoremtxs = [], [], []

    # adj pvalue of each gene for each group
    pvalmtx = pd.DataFrame( index= genes)
    # logfold change of each gene for each group
    logfoldmtx =  pd.DataFrame(index=genes)
    # zscore of each gene for each group
    scoremtx =  pd.DataFrame(index = genes)

    print(" Computing comparison: " + exp.obs['status'].unique())
    for ct in exp.obs['cell_type'].unique():
        print(ct)
        subset = exp[(exp.obs['cell_type']== ct)]

        # Check if there is a minimum population frquency of cell types - status to see if the analysis is viable,
        # otherwise skip the analysis and warn the user about what was skipped.
        counts = Counter(zip(subset.obs['cell_type'], subset.obs['status']))
        #minimum cell type frequency
        min_ctfreq= 10
        frequent_ct = [key for key, count in counts.items() if count >= 10]
        if len(frequent_ct) != 2:
            print("The cell type coudn't be used in this comparison due to more than two status categories for each celltype or becasue there were less than 10 cell for a group comparison")
            print(counts)
            continue


        #Differential express analysis celltype_status1 vs celltype_status2
        # Rank genes to caracterize genes
        sc.tl.rank_genes_groups(subset, groupby='status', reference=base_state, n_genes=subset.shape[1], key_added = '{}_DE'.format(ct), method='wilcoxon')

        #Write raw results of DE to file 
        genes = list(set(exp.var_names))
        gene2idx = {gene:i for i, gene in enumerate(genes)}
        pvalmtxs, logfoldmtxs, scoremtxs = [], [], []

        # loop through and fill up the matrix with pvalue, logfold and score
        iters =  zip(subset.uns['{}_DE'.format(ct)]['names'], 
                    subset.uns['{}_DE'.format(ct)]['pvals_adj'], 
                    subset.uns['{}_DE'.format(ct)]['logfoldchanges'], 
                    subset.uns['{}_DE'.format(ct)]['scores'])
        pvalrow = []
        logfoldrow=[]
        scorerow = []
        for gene, pval, logfold, score in iters:
            if gene[0] in list(gene2idx.keys()):
                pvalrow.append(pval[0])
                logfoldrow.append(logfold[0])
                scorerow.append(score[0])
            else :
                pvalrow.append(math.nan)
                logfoldrow.append(math.nan)
                scorerow.append(math.nan)
        res_name = diseaselab + "vs" + base_state +"_"+ct
        pvalmtx[res_name] = pvalrow
        logfoldmtx[res_name] = logfoldrow
        scoremtx[res_name] = scorerow

    # write matrix to file
    # This matrix contains a cloumn per cell type and rows per genes
    # It stores the results of DE of one celltype vs all other celltypes
    if os.path.isdir('./status_program') == False:
        os.mkdir("./status_program")
    save_dir="./status_program"
    pvalmtx.to_csv("%s/%s%svs%s_adjpval.csv"%(save_dir, data_name, diseaselab, base_state))
    logfoldmtx.to_csv("%s/%s%svs%s_logfold.csv"%(save_dir, data_name, diseaselab, base_state))
    scoremtx.to_csv("%s/%s%svs%s_zscore.csv"%(save_dir, data_name, diseaselab, base_state))

    ### Z - Score transformation ###
    # Transform zscore to x=-2log(zscore)
    scoremtx2 = scoremtx.clip(lower=0)
    scoremtx2 = pd.DataFrame(scipy.stats.norm.sf(scoremtx2), index=scoremtx2.index, columns=scoremtx2.columns)
    scoremtx2 = scoremtx2+1e-08
    scoremtx2 = -2*np.log(scoremtx2)
    # Minimization
    scoremtx2 = (scoremtx2 - scoremtx2.min())/(scoremtx2.max()-scoremtx2.min())
    # Save transfromed scores
    scoremtx2.to_csv("%s/%s_transgenescores.csv"%(save_dir, data_name))



    # Simple transformation
    # Transform zscore to x=-2log(zscore)
    scoremtxs = pd.DataFrame(scipy.stats.norm.sf(scoremtxs2), index=scoremtxs2.index, columns=scoremtxs2.columns)
    scoremtxs = scoremtxs2+1e-08
    scoremtxs = -2*np.log(scoremtxs)

    # Minimization
    scoremtxs = (scoremtxs - scoremtxs.min())/(scoremtxs.max()-scoremtxs.min())

    #Save transformed scores
    scoremtxs.to_csv("%s/%s_simpletransgenescores.csv"%(save_dir, data_name))