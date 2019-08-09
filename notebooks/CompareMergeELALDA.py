#!/usr/bin/env python
# coding: utf-8

# In[11]:


import sklearn
import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt
import seaborn as sns
import scanpy as sc
import anndata
import re
import importlib.util
spec = importlib.util.spec_from_file_location("ScanpyUtilsMT", os.path.expanduser("~/code/pollye/MTsc/utils/ScanpyUtilsMT.py"))
sc_utils = importlib.util.module_from_spec(spec)
spec.loader.exec_module(sc_utils)

subtractive=False
headpath='/wynton/scratch/mtschmitz/macaqueseq/'
adatapaths=[os.path.join(headpath,x) for x in  os.listdir(headpath)]
samplenames=[re.sub('_Out','',x) for x in  os.listdir(headpath)]

adatapaths=[x for x in adatapaths if 'PEC_YALE' in x]
samplenames=[x for x in samplenames if 'PEC_YALE' in x]
region=[x.split('_')[4] for x in samplenames]
from collections import defaultdict
regdict=defaultdict(list)
for r in list(set(region)):
    for i in range(len(region)):
        if region[i] == r:
            regdict[r].append(adatapaths[i])
print(regdict,flush=True)

for reg,item in regdict.items():
    sc.settings.figdir=os.path.expanduser('~/figs/'+str(subtractive)+'/merging'+reg)
    if not os.path.exists(sc.settings.figdir):
            os.makedirs(sc.settings.figdir)
    if not os.path.exists(os.path.join(regdict[reg][0],'outs/ambientsubtracted'+str(subtractive)+'.h5ad')):
        continue
    if not os.path.exists(os.path.join(regdict[reg][1],'outs/ambientsubtracted'+str(subtractive)+'.h5ad')):
        continue

    adatas=[]
    for p in regdict[reg]:
        a=sc.read_10x_mtx(os.path.join(p,'outs/filtered_feature_bc_matrix'),cache=True)
        sc.pp.filter_genes(a, min_cells=10,inplace=True)
        sc.pp.filter_cells(a,min_counts=5,inplace=True)
        adatas.append(a)
    adata=sc.AnnData.concatenate(*adatas)    
    
    
    adatas=[]
    for p in regdict[reg]:
        a=sc.read_h5ad(os.path.join(p,'outs/ambientsubtracted'+str(subtractive)+'.h5ad'))
        sc.pp.filter_genes(a, min_cells=10,inplace=True)
        sc.pp.filter_cells(a,min_counts=5,inplace=True)
        adatas.append(a)
    
    ambinotadata=sc.AnnData.concatenate(*adatas)

    cells=list(set(adata.obs.index) & set(ambinotadata.obs.index)) 
    genes=list(set(adata.var.index) & set(ambinotadata.var.index)) 
    
    ambinotadata._inplace_subset_obs(cells)
    ambinotadata._inplace_subset_var(genes)
    adata._inplace_subset_obs(cells)
    adata._inplace_subset_var(genes)
    sns.distplot((adata.X-ambinotadata.X).sum(1).A1,kde=False)
    plt.savefig(os.path.join(sc.settings.figdir,'CountDifferences.png'))
    plt.close()
    print(ambinotadata.var)
    b0=adata[adata.obs['batch']=='0',:].X.sum(0)/adata[adata.obs['batch']=='0',:].X.sum()
    b1=adata[adata.obs['batch']=='1',:].X.sum(0)/adata[adata.obs['batch']=='0',:].X.sum()
    print(np.corrcoef((b0-b1).A1,ambinotadata.var['lda_9-1']),flush=True)
    print(np.corrcoef((adata.X-ambinotadata.X).sum(0).A1,ambinotadata.var['lda_9-1']),flush=True)
    print(ambinotadata.var.iloc[(b0-b1).A1.argsort(),:])
    print(ambinotadata.var.iloc[(b1-b0).A1.argsort(),:])

    condit='Original'
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    #vg=sc.pp.highly_variable_genes(adata,adata.shape[0],inplace=False)
    #adata._inplace_subset_var(vg['highly_variable'])
    sc.pp.scale(adata, max_value=10)
    sc.pp.pca(adata)
    sc.pp.neighbors(adata)
    sc.tl.umap(adata)
    sc.tl.leiden(adata)
    sc.pl.umap(adata, color=['leiden','batch'],save=condit+"Leiden")
    sc_utils.marker_analysis(adata,variables=['leiden','batch'],markerpath=os.path.expanduser('~/markers.txt'),prefix=condit)
    condit='Ambiaint'
    sc.pp.normalize_total(ambinotadata, target_sum=1e4)
    sc.pp.log1p(ambinotadata)
    #ambinotadata._inplace_subset_var(vg['highly_variable'])
    sc.pp.scale(ambinotadata, max_value=10)
    sc.pp.pca(ambinotadata)
    sc.pp.neighbors(ambinotadata)
    sc.tl.umap(ambinotadata)
    sc.tl.leiden(ambinotadata)
    sc.pl.umap(ambinotadata, color=['leiden','batch'],save=condit+"Leiden")
    sc_utils.marker_analysis(ambinotadata,variables=['leiden','batch'],markerpath=os.path.expanduser('~/markers.txt'),prefix=condit) 
    pd.DataFrame([sklearn.metrics.silhouette_score(adata.X,adata.obs['leiden']) ,sklearn.metrics.silhouette_score(ambinotadata.X,adata.obs['leiden'])],index=['Original','Ambiaint']).to_csv(os.path.join(sc.settings.figdir,"SpectralScores.csv"))
    arand=sklearn.metrics.adjusted_rand_score(adata.obs['batch'],adata.obs['leiden'])
    aintrand=sklearn.metrics.adjusted_rand_score(ambinotadata.obs['batch'],ambinotadata.obs['leiden'])
    pd.DataFrame([arand,aintrand],index=['Original','Ambiaint']).to_csv(os.path.join(sc.settings.figdir,"BatchXLeidenRand.csv"))
    
    
