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

headpath='/wynton/scratch/mtschmitz/macaqueseq/'
adatapaths=[os.path.join(headpath,x) for x in  os.listdir(headpath)]
samplenames=[re.sub('_Out','',x) for x in  os.listdir(headpath)]
subtractive=True

for adatapath,samplename in zip(np.array(adatapaths),np.array(samplenames)):
    sc.settings.figdir=os.path.expanduser('~/figs/'+str(subtractive)+'/'+samplename)
    if not os.path.exists(sc.settings.figdir):
            os.makedirs(sc.settings.figdir)
    if not os.path.exists(os.path.join(adatapath,'outs/ambientsubtracted'+str(subtractive)+'.h5ad')):
        continue
    adata=sc.read_10x_mtx(os.path.join(adatapath,'outs/filtered_feature_bc_matrix'),cache=True)
    sc.pp.filter_genes(adata, min_cells=10,inplace=True)
    sc.pp.filter_cells(adata,min_counts=5,inplace=True)
    ambinotadata=sc.read_h5ad(os.path.join(adatapath,'outs/ambientsubtracted'+str(subtractive)+'.h5ad'))
    sc.pp.filter_genes(ambinotadata, min_cells=10,inplace=True)
    sc.pp.filter_cells(ambinotadata,min_counts=5,inplace=True)
    cells=list(set(adata.obs.index) & set(ambinotadata.obs.index)) 
    genes=list(set(adata.var.index) & set(ambinotadata.var.index)) 
    
    ambinotadata._inplace_subset_obs(cells)
    ambinotadata._inplace_subset_var(genes)
    adata._inplace_subset_obs(cells)
    adata._inplace_subset_var(genes)
    sns.distplot(np.log10((adata.X-ambinotadata.X).sum(1).A1+1),kde=False)
    plt.savefig(os.path.join(sc.settings.figdir,'CountDifferences.png'))
    plt.close()
    
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
    sc.pl.umap(adata, color=['leiden'],save=condit+"Leiden")
    sc_utils.marker_analysis(adata,variables=['leiden'],markerpath=os.path.expanduser('~/markers.txt'),subclass=True,prefix=condit)
    condit='Ambiaint'
    sc.pp.normalize_total(ambinotadata, target_sum=1e4)
    sc.pp.log1p(ambinotadata)
    #ambinotadata._inplace_subset_var(vg['highly_variable'])
    sc.pp.scale(ambinotadata, max_value=10)
    sc.pp.pca(ambinotadata)
    sc.pp.neighbors(ambinotadata)
    sc.tl.umap(ambinotadata)
    sc.tl.leiden(ambinotadata)
    sc.pl.umap(ambinotadata, color=['leiden'],save=condit+"Leiden")
    sc_utils.marker_analysis(ambinotadata,variables=['leiden'],markerpath=os.path.expanduser('~/markers.txt'),subclass=False,prefix=condit) 
    pd.DataFrame([sc_utils.get_cluster_metrics(adata,rands=[],silhouettes=['X_umap','X_pca']),sc_utils.get_cluster_metrics(ambinotadata,rands=[],silhouettes=['X_umap','X_pca'])],index=['Original','Ambiaint']).to_csv(os.path.join(sc.settings.figdir,"SilhouetteScores.csv"))
    
    
