# Opening the pickled results of NMF

# Load the object class
import scanpy as sc
from collections import Counter
import numpy as np
import pandas as pd 
import anndata
import os
import seaborn as sns
import scipy
from scipy import sparse
from time import time
from sklearn.utils import check_random_state
from sklearn.utils.extmath import safe_sparse_dot
from sklearn.decomposition import NMF
from matplotlib import pyplot as plt
import random
from scipy.sparse import csr_matrix
import pickle



_largenumber = 1E100
_smallnumber = 1E-10
_smallnumber2 = 1E-04

class JointNMF:
    def __init__(self, exp_mat_dir, status_colname, hlabel, dlabel , Wh=None, Hh=None, Wshh=None, Hshh=None, Wshd=None, Hshd=None, Wd=None, Hd=None, 
                 nh_components=10, nsh_components=10, nd_components=10, gamma=1, mu=None, numstarts=5):
        
        self.exp = sc.read_h5ad(exp_mat_dir)
        self.status_colname = status_colname
        self.exp = self.exp[(self.exp.obs[self.status_colname] == hlabel) | (self.exp.obs[self.status_colname] == dlabel)]
        scaleX = self.exp.X.min()

        self.Xh = self.exp[self.exp.obs[self.status_colname]==hlabel]
        self.Xd = self.exp[self.exp.obs[self.status_colname]==dlabel]

        # Trimmdown the matrixes to be the same cellcount
        # In later versions modify this to intead of picking cells randomply, pick cells
        # randomply in proportion to thier population count (in terma of data set, status and celltype)
        min_rows = min(self.Xh.shape[0], self.Xd.shape[0])
        if self.Xh.shape[0] > min_rows:
            # Truncate Wh and substitute the last row with its mean
            rows_to_delete = random.sample(range(self.Xh.shape[0]), self.Xh.shape[0]- min_rows)
            rows_to_keep = np.setdiff1d(np.arange(self.Xh.shape[0]), rows_to_delete)
            Xh_truncated = self.Xh[rows_to_keep]
            self.Xh = Xh_truncated
    
        if self.Xd.shape[0] > min_rows:
            # Truncate Wd and substitute the last row with its mean
            rows_to_delete = random.sample(range(self.Xd.shape[0]), self.Xd.shape[0]- min_rows)
            rows_to_keep = np.setdiff1d(np.arange(self.Xd.shape[0]), rows_to_delete)
            Xd_truncated = self.Xd[rows_to_keep]
            self.Xd = Xd_truncated

        #######
        # When gnerating Latent factors data, first give some information
        # on the data that is being processed
        print("Generating celltype programs for the file:" + exp_mat_dir)
        self.data_name = exp_mat_dir[exp_mat_dir.rfind("/") + 1:(exp_mat_dir.rfind(".h5ad"))]
        print("Automatically chosen data name:  " + self.data_name)
        print("\nDescription of the data:")
        print("Healthy expression matrix : " + str(self.Xh.shape[0]) + " cells | " + str(self.Xh.shape[1]) + " genes")
        print("Disease expression matrix : " +str(self.Xd.shape[0]) + " cells | " + str(self.Xd.shape[1]) + " genes")
        print("cell state groups: " + str(list(self.exp.obs[status_colname].unique())))
        print("Generating joint NMF for " + hlabel + "and " + dlabel)
        ######

        Xh = self.Xh.X
        Xd = self.Xd.X
        if scaleX < 0 :
            Xh = Xh + scaleX
            Xd = Xd + scaleX
        
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
        
        # Healthy programs
        # generate helathy - latent factor contribution
        if Wh is None:
            self.Wh = scipy.sparse.csr_matrix(self.nmfh.transform(self.Xh))
        else:
            if (Wh.shape != (self.Xh.shape[0], self.nh_components + self.nsh_components)):
                raise ValueError("Initial Wh has wrong shape.")
            self.Wh = np.copy(Wh)
        # generate healthy - latent factor contribution
        if Hh is None:
            self.Hh = scipy.sparse.csr_matrix(self.nmfh.components_)
        else:
            if (Hh.shape != (self.nh_components + self.nsh_components, self.Xh.shape[1])):
                raise ValueError("Initial Wh has wrong shape.")
            self.Hh = np.copy(Hh) 
        
        # Disease programs
        # generate disease - latent factor contribution
        if Wd is None:
            self.Wd = scipy.sparse.csr_matrix(self.nmfd.transform(self.Xd))
        else:
            if (Wd.shape != (self.Xd.shape[0], self.nd_components + self.nsh_components)):
                raise ValueError("Initial Wd has wrong shape.")
            self.Wd = np.copy(Wd)
        # generate disease - latent factor contrbution
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
        
        # Computing reconstruction error (How simmilar are the matrixes to the original ones)
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


        #Save data

        write_mat_jointNMF(self.Wd, "Wd", self.data_name)
        write_mat_jointNMF(self.Hd, "Hd", self.data_name)

        write_mat_jointNMF(self.Wh, "Wh", self.data_name)
        write_mat_jointNMF(self.Hh, "Hh", self.data_name)

        write_mat_jointNMF(self.Wshd, "Wshd", self.data_name)
        write_mat_jointNMF(self.Wshd, "Wshh", self.data_name)

        print( "(final chi , time used)")

        return (chi2, time_used)
    

def write_mat_jointNMF(mtx, name, data_name):
    # Write and save results in folder
    save_dir="./process_program" + data_name
    if os.path.isdir(save_dir) == False:
        os.mkdir(save_dir)
    df = mtx.toarray()
    df = pd.DataFrame(df, columns = ["NMF" + str(i) for i in range(df.shape[1])])
    df.to_csv("%s/%s_%s.csv"%(save_dir,name))




def get_file_name(file_dir):
    last_slash_index = file_dir.rfind('/')
    pkl_index = file_dir.rfind('.pkl')
    substring = file_dir[last_slash_index + 1 : pkl_index]
    return substring

##########################################
    # SET UP END
##########################################




##########################################
file_dir= "ADvsCrl_pp.pkl"
case= 'AD'
file_name = get_file_name(file_dir)
with open(file_dir, 'rb') as picklefile:
    this_NMF = pickle.load(picklefile)

# Create initial correlation matrix
correlations = []
for i in range(this_NMF.Wh.shape[1]):
    correlation = []
    for j in range(this_NMF.Wd.shape[1]):
        corr = np.corrcoef(this_NMF.Wh[:,i].T.todense(), this_NMF.Wd[:,j].T.todense())
        correlation.append(corr[1][0])
    correlations.append(correlation)
correlations = np.array(correlations)
plt.imshow(correlations)
plt.xlabel("Disease Programs")
plt.ylabel("Healthy Programs")
plt.title("Correlation Heatmap "+ case + " vs Control")
plt.colorbar()
plt.savefig("correlation_initial_heatmap"+ file_name +".png")
plt.show()




file_dir= "COVIDvsCrl_pp.pkl"
case= 'COVID'
file_name = get_file_name(file_dir)
with open(file_dir, 'rb') as picklefile:
    this_NMF = pickle.load(picklefile)

# Create initial correlation matrix
correlations = []
for i in range(this_NMF.Wh.shape[1]):
    correlation = []
    for j in range(this_NMF.Wd.shape[1]):
        corr = np.corrcoef(this_NMF.Wh[:,i].T.todense(), this_NMF.Wd[:,j].T.todense())
        correlation.append(corr[1][0])
    correlations.append(correlation)
correlations = np.array(correlations)
plt.imshow(correlations)
plt.xlabel("Disease Programs")
plt.ylabel("Healthy Programs")
plt.title("Correlation Heatmap "+ case + " vs Control")
plt.colorbar()
plt.savefig("correlation_initial_heatmap"+ file_name +".png")
plt.show()



file_dir= "CasevsCtrl_pp.pkl"
case= 'Case'
file_name = get_file_name(file_dir)
with open(file_dir, 'rb') as picklefile:
    this_NMF = pickle.load(picklefile)

# Create initial correlation matrix
correlations = []
for i in range(this_NMF.Wh.shape[1]):
    correlation = []
    for j in range(this_NMF.Wd.shape[1]):
        corr = np.corrcoef(this_NMF.Wh[:,i].T.todense(), this_NMF.Wd[:,j].T.todense())
        correlation.append(corr[1][0])
    correlations.append(correlation)
correlations = np.array(correlations)
plt.imshow(correlations)
plt.xlabel("Disease Programs")
plt.ylabel("Healthy Programs")
plt.title("Correlation Heatmap "+ case + " vs Control")
plt.colorbar()
plt.savefig("correlation_initial_heatmap"+ file_name +".png")
plt.show()