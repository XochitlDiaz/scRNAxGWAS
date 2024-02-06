# 1 Feb 2024
# Simulated data for scRNAseq data
# XochitlDiaz

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
#%%
### Counts Matrix

genes = ["CD11B", "CD14", "CD16", "CD40", "CD45", "CD80", "CD68", "CD115", 
         "CX3CR1", "TMEM119", "FCER1G", "FCRLS", "SLC2A5", "P2Y12", "TSPO", 
         "ARG1", "FIZZ1", "MMP9", "MMP12", "APOE4"]


cells = ['AAA', 'AAC', 'AAT', 'AAG', 'ACA', 'ACC', 'ACT', 'ACG', 'ATA', 'ATC', 'ATT', 'ATG', 'AGA', 'AGC', 'AGT', 'AGG', 'CAA', 'CAC', 'CAT', 'CAG', 'CCA', 'CCC', 'CCT', 'CCG', 'CTA', 'CTC', 'CTT', 'CTG', 'CGA', 'CGC', 'CGT', 'CGG', 'TAA', 'TAC', 'TAT', 'TAG', 'TCA', 'TCC', 'TCT', 'TCG', 'TTA', 'TTC', 'TTT', 'TTG', 'TGA', 'TGC', 'TGT', 'TGG', 'GAA', 'GAC', 'GAT', 'GAG', 'GCA', 'GCC', 'GCT', 'GCG', 'GTA', 'GTC', 'GTT', 'GTG', 'GGA', 'GGC', 'GGT', 'GGG']

prefix = "AAACAT"

suffix = "-1"

cells = [prefix + element + suffix for element in cells]

df = pd.DataFrame(0,index = cells, columns= genes)

# Defined random values for rows 0:43
for i in range(0,len(df.index)-20):
    df.iloc[[i]] = np.random.normal(loc=5, scale=2, size= len(df.columns))

# Downreg case 1
df.iloc[63-5:64, 0:6] = -5
# Upreg case 1
df.iloc[63-5:64, 5:11]= 5

# Up reg case 2
df.iloc[57-5:58, 0:6] = 5
# Down reg case 3
df.iloc[57-5:58, 5:11]= -5

df.to_csv("./dummy_data/counts.csv", index= False, header= False)

# %%
#Index names
df.index
# single cell value
df.iloc[[0],[7]]
#Row values
df.iloc[[0]]
# row values as a list 
df.iloc[0]
# values of two rows
df.iloc[[0,4]]
# values of a range of rows
df.iloc[63-5:63]
# all valuense of column 0
df.iloc[:,0]

# Empty dataframe
df2 = pd.DataFrame()
# Define data as list of lists
data = [['tom', 10], ['nick', 15], ['juli', 14]]
# Create the pandas DataFrame
df2 = pd.DataFrame(data, columns=['Name', 'Age'])

# initialize data of lists.
data = {'Name': ['Tom', 'nick', 'krish', 'jack'],
        'Age': [20, 21, 19, 18]} 
# Create DataFrame
df2 = pd.DataFrame(data)

#Write DataFrame to csv
df.to_csv("./dummy_data/counts.csv")

#%%
### Meta Data Matrix
age = np.random.normal(loc=65, scale= 5, size=len(df.index))
md_data = {"Patient": np.repeat(["P1","P2","P3","P4","P5","P6","P7"],[10,10,10,10,12,6,6]),
            "Status": np.repeat(["HC1", "HC2", "Case1", "Case2"],[44,8,6,6]),
            "Age": np.around(age,0),
            "Celltype": np.repeat("mic",len(df.index)),  }
md= pd.DataFrame(md_data,index= cells)

md.to_csv("./dummy_data/metadata.csv")
