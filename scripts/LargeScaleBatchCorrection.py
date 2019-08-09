import sklearn
import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt
import seaborn as sns
import scanpy as sc
import anndata
import re
import bbknn
import importlib.util
spec = importlib.util.spec_from_file_location("ScanpyUtilsMT", os.path.expanduser("~/code/pollye/MTsc/utils/ScanpyUtilsMT.py"))
sc_utils = importlib.util.module_from_spec(spec)
spec.loader.exec_module(sc_utils)
sc.settings.figdir=os.path.expanduser('~/figs/LargeScaleBatchCorrect')
headpath='/wynton/scratch/mtschmitz/macaqueseq/'
adatapaths=[os.path.join(headpath,x) for x in  os.listdir(headpath)]
samplenames=[re.sub('_Out','',x) for x in  os.listdir(headpath)]
adatapaths=[os.path.join(x,'outs/filtered_feature_bc_matrix') for x in adatapaths]
adatas=[]
for adatapath,samplename in zip(np.array(adatapaths),np.array(samplenames)):
    print(samplename,flush=True)
    a=sc.read_10x_mtx(adatapath)
    cells=np.random.choice(a.obs.index,min(a.shape[0],500),replace=False)
    a._inplace_subset_obs(cells)
    adatas.append(a)

adata=sc.AnnData.concatenate(*adatas)
adata.var_names_make_unique()

mito_genes = [name for name in adata.var_names if name in ['ND1','ND2','ND4L','ND4','ND5','ND6','ATP6','ATP8','CYTB','COX1','COX2','COX3'] or 'MT-' in name]
ribo_genes = [name for name in adata.var_names if name.startswith('RPS') or name.startswith('RPL') ]

adata.obs['percent_mito'] = np.sum(
    adata[:, mito_genes].X, axis=1).A1 / np.sum(adata.X, axis=1).A1
adata.obs['percent_ribo'] = np.sum(
    adata[:, ribo_genes].X, axis=1).A1 / np.sum(adata.X, axis=1).A1
# add the total counts per cell as observations-annotation to adata
adata.obs['n_counts'] = adata.X.sum(axis=1).A1
#scanpy.api.pp.filter_genes_dispersion(adata,n_top_genes=np.sum(np.sum(adata.X, axis=0)>0))
sc.pp.filter_genes(adata, min_cells=20,inplace=True)
sc.pp.filter_cells(adata,min_genes=400)
#sc.pl.violin(adata, ['n_genes', 'n_counts', 'percent_mito'],jitter=0.4,multi_panel=True)
print(adata,flush=True)
adataraw=adata.copy()
sc.pp.normalize_per_cell(adata, counts_per_cell_after=1e4)
sc.pp.log1p(adata)
adata.raw = adata
normalizedadata=adata.copy()
sc.pp.highly_variable_genes(adata, n_top_genes=5000)
sc.pp.scale(adata, max_value=10)
sc.pp.pca(adata)
adatabbknn=adata.copy()
sc.pp.neighbors(adata)
sc.tl.umap(adata)
sc.tl.leiden(adata)
sc.pl.umap(adata,color=['leiden','batch'],save='_nocorrect')
print('Metrics',flush=True)
df=pd.DataFrame([sc_utils.get_cluster_metrics(adata,rands=['batch'])])
print('RobotAlex',flush=True)
sc_utils.marker_analysis(adata,variables=['leiden'],markerpath=os.path.expanduser('~/markers.txt'),subclass=True,prefix='nocorrect')
print('Markers',flush=True)
sc_utils.log_reg_diff_exp(adata,'nocorrect')
print('RegressingOut',flush=True)
adata.X = adata.raw.X
sc.pp.normalize_per_cell(adata, counts_per_cell_after=1e4)
sc.pp.log1p(adata)
sc.pp.regress_out(adata,'batch')
sc.pp.highly_variable_genes(adata, n_top_genes=5000)
sc.pp.scale(adata, max_value=10)
sc.pp.pca(adata)
sc.pp.neighbors(adata)
sc.tl.umap(adata)
sc.tl.leiden(adata)
df=df.append(sc_utils.get_cluster_metrics(adata,['batch']),ignore_index=True)
sc_utils.marker_analysis(adata,variables=['leiden'],markerpath=os.path.expanduser('~/markers.txt'),subclass=True,prefix='regressed')
sc_utils.log_reg_diff_exp(adata,'regressed')
sc.pl.umap(adata,color=['leiden','batch'],save='_regressed')
df.to_csv(os.path.join(sc.settings.figdir,"Metrics.csv"))
adata.X=adata.raw.X 
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, n_top_genes=5000)
sc.pp.combat(adata, key='batch')
sc.pp.scale(adata, max_value=10)
sc.pp.pca(adata)
sc.pp.neighbors(adata)
sc.tl.umap(adata)
sc.tl.leiden(adata)
df=df.append(sc_utils.get_cluster_metrics(adata,['batch']),ignore_index=True)
sc_utils.marker_analysis(adata,variables=['leiden'],markerpath=os.path.expanduser('~/markers.txt'),subclass=True,prefix='combat')
sc_utils.log_reg_diff_exp(adata,'combat')
bbknn.bbknn(adatabbknn)
sc.tl.umap(adata)
sc.tl.leiden(adata)
df=df.append(sc_utils.get_cluster_metrics(adata,['batch']),ignore_index=True)
sc_utils.marker_analysis(adata,variables=['leiden'],markerpath=os.path.expanduser('~/markers.txt'),subclass=True,prefix='bbknn')
sc_utils.log_reg_diff_exp(adata,'bbknn')
sc.pl.umap(adata,color=['leiden','batch'],save='_bbknn')
print(df)
df.to_csv(os.path.join(sc.settings.figdir,"Metrics.csv"))

sc.pp._dca.dca(adataraw)
sc.pp.neighbors(adataraw,rep='X_dca')
sc.tl.leiden(adataraw)
df=df.append(sc_utils.get_cluster_metrics(adataraw,['batch'],silhouettes=['X_pca','X_dca']),ignore_index=True)
sc_utils.marker_analysis(adataraw,variables=['leiden'],markerpath=os.path.expanduser('~/markers.txt'),subclass=True,prefix='dca')
sc_utils.log_reg_diff_exp(adataraw,'dca')
sc.pl.umap(adata,color=['leiden','batch'],save='_dca')

sc.pp.pca(adataraw)
sc.pp.neighbors(adataraw)
sc.tl.leiden(adataraw)
df=df.append(sc_utils.get_cluster_metrics(adataraw,['batch']),ignore_index=True)
sc_utils.marker_analysis(adataraw,variables=['leiden'],markerpath=os.path.expanduser('~/markers.txt'),subclass=True,prefix='dcapca')
sc_utils.log_reg_diff_exp(adataraw,'dcapca')
sc.pl.umap(adata,color=['leiden','batch'],save='_dcapca')

print(df)
df.to_csv(os.path.join(sc.settings.figdir,"Metrics.csv"))