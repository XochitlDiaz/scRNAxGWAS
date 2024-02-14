# Generate Cell type Programs
# 9 Feb 2024
# XochitlDiaz

# This function will generate Cell Type programs out 
# of a previously QCed .h5ad file

import scanpy as sc
from collections import Counter
import numpy as np
import pandas as pd 
import anndata
import os
import seaborn as sns
from scipy import sparse
from time import time
from sklearn.utils import check_random_state
from sklearn.utils.extmath import safe_sparse_dot
from matplotlib import pyplot as plt


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


########## TEST ###############
scdatasets = {'simdata': ("./dummy_data/simdata.h5ad","Patient") }

def genCellularProcessPrograms(exp_mat_dir, ct_colname, status_colname, sample_colname):
    #read in data
    exp = sc.read_h5ad(exp_mat_dir) # tissuedata=exp
    #change referenced column names for anndata manipulation
    exp.obs.rename(columns={ct_colname:'cell_type', status_colname:'status', sample_colname:'patient_id'},inplace=True)
    
    num_celltypes = len(exp.obs["cell_type"].unique())
    #######
    # When gnerating Celltype Data, first give some information
    # on the data that is being processed
    print("Generating celltype programs for the file:" + exp_mat_dir)
    data_name = exp_mat_dir[exp_mat_dir.rfind("/") + 1:(exp_mat_dir.rfind(".h5ad"))] # data_name = tissue
    print("Automatically chosen data name:  " + data_name)
    print("\nDescription of the data:")
    print( str(num_celltypes) + " cells | " + str(exp.shape[1]) + " genes | "+ str(len(set(exp.obs['patient_id']))) + " samples")
    print(str(len(set(exp.obs['cell_type']))) + " cell types: " + str(list(exp.obs['cell_type'].unique())))
    print("cell state groups: " + str(list(exp.obs['status'].unique())))
    print("\nComputing " + str(num_celltypes + 10) + " cellular processes")
    ######

    model = NMF(n_components=num_celltypes+10, init='random', random_state=0)
    sc.pp._highly_variable_genes(exp, n_top_genes=2000)
    exp = exp[:,exp.var['highly_variable']]
    exp.layers['counts'] = exp.X
    X = exp.layers['counts']
    # NMF only accepts non negative values
    if X.min()< 0:
        X = X + min(X)
    X = X/np.max(X)
    W = model.fit_transform(X)
    W = pd.DataFrame(W, index = exp.obs_names, columns = ['NMF_%d'%i for i in range(W.shape[1])])
    H = pd.DataFrame(model.components_.T, index=exp.var_names, columns=['NMF_%d'%i for i in range(model.components_.shape[0])])

    if not os.path.exists("./process_program"):
        os.mkdir("./process_program")

    save_dir = "./process_program"
    W.to_csv("%s/%s_cellprograms.csv"%(save_dir,data_name))
    H.to_csv("%s/%s_geneprograms.csv"%(save_dir,data_name))

    for i in range(W.shape[1]):
        exp.obs['NMF_%d'%i] = W['NMF_%d'%i]
    
    if 'umap' not in exp.obsm:
        sc.pp.neighbors(exp)
        sc.tl.umap(exp)
    sc.pl.umap(exp, color=['NMF_%d'%i for i in range(W.shape[1])], color_map='Reds', save="%s/%s_healthyprograms.pdf"%(save_dir,data_name))
    sc.pl.umap(exp, color=['cell_type'], ncols=1, save="%s/%s_healthyUMAP.pdf"%(save_dir,data_name))
    
    
    # same two plots as above but they are joint together on the same margin
    sc.pl.umap(exp, color=['NMF_%d'%i for i in range(W.shape[1])]+['cell_type'], 
                   color_map='Reds', save='%s/%s_healthy.png'%(save_dir,data_name))
    

    #Compute correlation between matrices x and y 
    # Where x is the expression of cells and y is the contribution of latent factors to a cell
    x = scipy.sparse.csc_matrix(exp.X)
    #x = x[:,0:-6]
    nrow,ncol = x.shape
    correlations = []
    for i in range(W.shape[1]):
        y = W['NMF_%d'%i].values
        yy = y - y.mean()
        xm = x.mean(axis=0).A.ravel()
        ys = yy / np.sqrt(np.dot(yy, yy))
        xs = np.sqrt(np.add.reduceat(x.data**2, x.indptr[:-1]) - nrow*xm*xm)

        correl2 = np.add.reduceat(x.data * ys[x.indices], x.indptr[:-1]) / xs
        correlations.append(correl2)
    corr = pd.DataFrame(np.vstack(correlations).T, columns=W.columns, index=exp.var_names)
    corr.to_csv('%s/%s_correlation_cellprograms.csv'%(save_dir,data_name))


    # Compress latent factors that had a correlation >0.8
    W0 = W
    W = W/W.quantile(0.99, axis=0)
    for i in range(W.shape[1]):
        exp.obs['NMF_%d'%i] = W['NMF_%d'%i]

    topicdata = sc.AnnData(W.values)
    topicdata.obs_names = W.index
    topicdata.var_names = W.columns
    topicdata.obs['cell_type'] = exp.obs['cell_type']

    sc.pl.matrixplot(topicdata, var_names=topicdata.var_names, groupby='cell_type', standard_scale='var' )
    
    W['cell_type'] = exp.obs['cell_type']

    W = W.groupby('cell_type').agg('mean')
    
    nmf_specificity = (W > .3).sum(0)
    display(W)
    for i in range(W.shape[0]):
        for j in range(W.shape[1]):
            if W.values[i,j] > .65 and nmf_specificity[j]==1:
                print(W.index[i], W.columns[j])
    
    
    #H is the gene programs maatrix
    print(H.shape)
    for i in range(H.shape[1]):
        print(i)
        print(",".join(H.sort_values(by='NMF_%d'%i, ascending=False).index[0:50]))
        print(",".join(H.sort_values(by='NMF_%d'%i, ascending=False).index[0:50]))
        
        print(len(set(H.sort_values(by='NMF_%d'%i, ascending=False).index[0:500]).intersection(H.sort_values(by='NMF_%d'%i, ascending=False).index[0:500])))
 
###############################


# Joint NMF
        
_largenumber = 1E100
_smallnumber = 1E-10
_smallnumber2 = 1E-04

class JointNMF:
    def __init__(self, Xh, Xd, Wh=None, Hh=None, Wshh=None, Hshh=None, Wshd=None, Hshd=None, Wd=None, Hd=None, 
                 nh_components=10, nsh_components=10, nd_components=10, gamma=1, mu=None, numstarts=5):
        self.Xh = scipy.sparse.csr_matrix(Xh).copy()
        self.Xd = scipy.sparse.csr_matrix(Xd).copy()
        self.Xh = self.Xh/np.max(self.Xh)
        self.Xd = self.Xd/np.max(self.Xd)

        # initializations
        self.nh_components = nh_components
        self.nsh_components = nsh_components
        self.nd_components = nd_components
        
        self.maxiters = 1000
        self.tol = _smallnumber
        
        self.gamma = gamma
        
        
        # initialize the matrix using result from standard NMF for healthy programs
        min_reconstruction_err = _largenumber
        bestnmfh = None
        for i in range(numstarts):
            nmfh = NMF(n_components = nsh_components + nh_components)
            nmfh.fit_transform(self.Xh)
            if min_reconstruction_err > nmfh.reconstruction_err_:
                min_reconstruction_err = nmfh.reconstruction_err_
                bestnmfh = nmfh
        self.nmfh = bestnmfh
        
        print("Initialized healthy with reconstruction error: ", min_reconstruction_err)
        
        # initialize the matrix using result from standard NMF for disease programs
        min_reconstruction_err = _largenumber
        bestnmfd = None
        for i in range(numstarts):
            nmfd = NMF(n_components = nsh_components + nd_components)
            nmfd.fit_transform(self.Xd)
            if min_reconstruction_err > nmfd.reconstruction_err_ :
                min_reconstruction_err = nmfd.reconstruction_err_
                bestnmfd = nmfd
        self.nmfd = bestnmfd
        
        print("Initialized disease with reconstruction error: ", min_reconstruction_err)
        
        # healthy programs
        if Wh is None:
            self.Wh = scipy.sparse.csr_matrix(self.nmfh.transform(self.Xh))
        else:
            if (Wh.shape != (self.Xh.shape[0], self.nh_components + self.nsh_components)):
                raise ValueError("Initial Wh has wrong shape.")
            self.Wh = np.copy(Wh)
            
        if Hh is None:
            self.Hh = scipy.sparse.csr_matrix(self.nmfh.components_)
        else:
            if (Hh.shape != (self.nh_components + self.nsh_components, self.Xh.shape[1])):
                raise ValueError("Initial Wh has wrong shape.")
            self.Hh = np.copy(Hh) 
        
        # disease programs
        if Wd is None:
            self.Wd = scipy.sparse.csr_matrix(self.nmfd.transform(self.Xd))
        else:
            if (Wd.shape != (self.Xd.shape[0], self.nd_components + self.nsh_components)):
                raise ValueError("Initial Wd has wrong shape.")
            self.Wd = np.copy(Wd)
        
        if Hd is None:
            self.Hd = scipy.sparse.csr_matrix(self.nmfd.components_)
        else:
            if (Hd.shape != (self.nd_components + self.nsh_components, self.Xd.shape[1])):
                raise ValueError("Initial Wd has wrong shape.")
            self.Hd = np.copy(Hd)

        reorder_healthy_idxs, reorder_disease_idxs = self.align_matrices()
        self.Wh_orig = self.Wh[:, reorder_healthy_idxs].copy()
        self.Wh = self.Wh[:, reorder_healthy_idxs]
        
        self.Hh_orig = self.Hh[reorder_healthy_idxs, :].copy()
        self.Hh = self.Hh[reorder_healthy_idxs, :]
        
        self.Wd_orig = self.Wd[:, reorder_disease_idxs].copy()
        self.Wd = self.Wd[:, reorder_disease_idxs]
        
        self.Hd_orig = self.Hd[reorder_disease_idxs, :].copy()
        self.Hd = self.Hd[reorder_disease_idxs, :]
        
        healthy_reconstruction = sparse.linalg.norm(self.Xh - self.Wh.dot(self.Hh), ord='fro')
        disease_reconstruction = sparse.linalg.norm(self.Xd - self.Wd.dot(self.Hd), ord='fro')
        print ("the reconstruction value for healthy and disease: %.5f, %.5f"%(healthy_reconstruction, disease_reconstruction))
        
        correlations = []
        for i in range(self.Wh.shape[1]):
            correlation = []
            for j in range(self.Wd.shape[1]):
                corr = np.corrcoef(self.Wh[:,i].T.todense(), self.Wd[:,j].T.todense())
                correlation.append(corr[1][0])
            correlations.append(correlation)
        correlations = np.array(correlations)
        plt.imshow(correlations)
        plt.xlabel("Disease Programs")
        plt.ylabel("Healthy Programs")
        plt.show()

        # option for user input mu or estimated mu 
        if mu:
            self.mu = mu
        else:
            healthy_diff = 0.5*sparse.linalg.norm(self.Xh - safe_sparse_dot(self.Wh, self.Hh), ord='fro')**2
            disease_diff = 0.5*sparse.linalg.norm(self.Xd - safe_sparse_dot(self.Wd, self.Hd), ord='fro')**2
            denominator = sparse.linalg.norm(self.Wh, ord='fro')**2 + sparse.linalg.norm(self.Wd, ord='fro')**2
            self.mu = (healthy_diff + disease_diff)/denominator
    
    def align_matrices(self):
        correlations = []
        for i in range(self.Wh.shape[1]):
            correlation = []
            for j in range(self.Wd.shape[1]):
                corr = np.corrcoef(self.Wh[:,i].T.todense(), self.Wd[:,j].T.todense())
                correlation.append(corr[1][0])
            correlations.append(correlation)
        correlations = np.array(correlations)
        correlations.shape

        reorder_healthy_idxs = []
        reorder_disease_idxs = []
        ct = 0

        while ct < min(correlations.shape[0], correlations.shape[1]):
            argmax = np.argmax([correlations[i, j] for i, j in enumerate(correlations.argmax(axis=1))])
            i = range(correlations.shape[0])[argmax]
            j = correlations.argmax(axis=1)[argmax]
            reorder_healthy_idxs.append(i)
            reorder_disease_idxs.append(j)
            correlations[i,:] = -2
            correlations[:,j] = -2
            ct = ct + 1

        reorder_healthy_idxs = reorder_healthy_idxs + list(set(range(correlations.shape[0])).difference(set(reorder_healthy_idxs)))
        reorder_disease_idxs = reorder_disease_idxs + list(set(range(correlations.shape[1])).difference(set(reorder_disease_idxs)))
        return reorder_healthy_idxs, reorder_disease_idxs
    
    @property
    def cost(self):
        self.Wh = self.Wh.tocsr()
        self.Wd = self.Wd.tocsr()
        self.Wshh = self.Wh[:,:self.nsh_components]
        self.Wshd = self.Wd[:,:self.nsh_components]
        
        diff1 = 0.5*sparse.linalg.norm(self.Xh - safe_sparse_dot(self.Wh, self.Hh), ord='fro')**2
        diff2 = 0.5*sparse.linalg.norm(self.Xd - safe_sparse_dot(self.Wd, self.Hd), ord='fro')**2
        diff3 = (self.mu/2)*(sparse.linalg.norm(self.Wh, ord='fro')**2) + (self.mu/2)*(sparse.linalg.norm(self.Wd, ord='fro')**2)
        diff4 = (self.gamma/2)*sparse.linalg.norm(self.Wshh-self.Wshd, ord='fro')**2
        chi2 = diff1 + diff2 + diff3 + diff4
        return chi2
    
    def solve(self, W_only=False, H_only=False, sparsemode=False, maxiters=None, tol=None):
        
        
        # initialize the state
        t0 = time()        
        niter = 0
        chi2 = self.cost
        oldchi2 = _largenumber
        maxiters = maxiters if maxiters else self.maxiters
        while (niter < maxiters) and np.abs((oldchi2-chi2)/oldchi2) > self.tol: #((oldchi2-chi2)/chi2 > self.tol):
            # update Wh
            scale2=np.append((self.gamma+self.mu)*np.ones(self.nsh_components), (self.mu)*np.ones(self.nh_components))
            self.Wshd = self.Wd[:,:self.nsh_components]
            Wh_up1 = safe_sparse_dot(self.Xh, self.Hh.T)
            Wshd_transform = self.Wshd.multiply(scipy.sparse.csr_matrix(self.gamma*np.ones((self.Xh.shape[0], self.nsh_components))))
            zeros = scipy.sparse.csr_matrix(np.zeros((self.Xh.shape[0], self.nh_components)))
            Wh_up2 = scipy.sparse.hstack((Wshd_transform, zeros)) #+ _smallnumber2
            Wh_down = safe_sparse_dot(self.Wh, safe_sparse_dot(self.Hh, self.Hh.T)) + safe_sparse_dot(self.Wh, np.diag(scale2)) #+ _smallnumber2
            self.Wh = self.Wh.multiply((Wh_up1 + Wh_up2)/Wh_down).tocsr()
            
            # update Hh
            Hh_up = safe_sparse_dot(self.Wh.T, self.Xh) #+ _smallnumber2
            Hh_down = safe_sparse_dot(safe_sparse_dot(self.Wh.T, self.Wh), self.Hh) #+ _smallnumber2
            self.Hh = self.Hh.multiply(Hh_up/Hh_down).tocsr()
            
            # update Wd
            scale2=np.append((self.gamma+self.mu)*np.ones(self.nsh_components), self.mu*np.ones(self.nd_components))
            self.Wshh = self.Wh[:,:self.nsh_components]
            Wd_up1 = safe_sparse_dot(self.Xd, self.Hd.T)
            Wshh_transform = self.Wshh.multiply(scipy.sparse.csr_matrix(self.gamma*np.ones((self.Xd.shape[0], self.nsh_components))))
            zeros = np.zeros((self.Xd.shape[0], self.nd_components))
            Wd_up2 = scipy.sparse.hstack((Wshh_transform, zeros)) #+ _smallnumber2
            Wd_down = safe_sparse_dot(self.Wd, safe_sparse_dot(self.Hd, self.Hd.T)) + safe_sparse_dot(self.Wd, np.diag(scale2)) #+ _smallnumber2 
            self.Wd = self.Wd.multiply((Wd_up1 + Wd_up2)/Wd_down).tocsr()
                        
            # update Hd
            Hd_up = safe_sparse_dot(self.Wd.T, self.Xd) #+ _smallnumber2
            Hd_down = safe_sparse_dot(safe_sparse_dot(self.Wd.T, self.Wd), self.Hd) #+ _smallnumber2
            self.Hd = self.Hd.multiply(Hd_up/Hd_down).tocsr()
                        
            # chi2
            oldchi2 = chi2
            chi2 = self.cost
            
            # Some quick check. May need its error class ...
            if (not np.isfinite(chi2)):
               raise ValueError("NMF construction failed, likely due to missing data")

            if (np.mod(niter, 10)==0):
                print("Current Chi2={0:.4f}, Previous Chi2={1:.4f}, \
                      Change={2:.4f}% @ niters={3}".format(chi2,oldchi2,((oldchi2-chi2)/oldchi2)*100.,niter), flush=True)
                """correlations = []
                for i in range(self.Wh.shape[1]):
                    correlation = []
                    for j in range(self.Wd.shape[1]):
                        corr = np.corrcoef(self.Wh[:,i].T.todense(), self.Wd[:,j].T.todense())
                        correlation.append(corr[1][0])
                    correlations.append(correlation)
                correlations = np.array(correlations)
                plt.imshow(correlations)
                plt.xlabel("Disease Programs")
                plt.ylabel("Healthy Programs")
                plt.show()"""

            niter += 1
            if (niter == self.maxiters):
                print("Iteration in re-initialization reaches maximum number = {0}".format(niter), flush=True)
                
        time_used = (time()-t0)/60.
        print("Took {0:.3f} minutes to reach current solution.".format(time_used), flush=True)

        return (chi2, time_used)
    







######## TEST ########
exp.obs['status'].unique()
scaleX = exp.X.min()
healthy_exp = exp[exp.obs['status']=='HC']
disese_exp = exp[exp.obs['status']=='Case']

Xh = healthy_exp.X
Xd = healthy_exp.X
if scaleX < 0 :
    Xh = Xh + scaleX
    Xd = Xd + scaleX
trail = JointNMF(Xh=Xh, Xd=Xd)

trail.align_matrices()
trail.cost
trail.solve()