import os
import scanpy
import anndata
import scanpy.api as sc
import pandas as pd
import numpy as np
import scipy
from scipy import stats
import re
import sklearn
import matplotlib
import matplotlib.pyplot as plt
matplotlib.use('Agg')
from umi_tools import umi_methods
from collections import Counter
import logging as logg
import random
import seaborn
import sys
import shutil

#Load my pipeline functions
import importlib.util
spec = importlib.util.spec_from_file_location("ScanpyUtilsMT", os.path.expanduser("~/ScanpyUtilsMT.py"))
sc_utils = importlib.util.module_from_spec(spec)
spec.loader.exec_module(sc_utils)
#Get lots of feedback from scanpy
logg.basicConfig(level=logg.INFO)
sc.settings.verbosity=5

refname="refdata-celranger-mmul8-toplevel"

sc.settings.autosave=True
sc.settings.autoshow=False

basepath=os.path.expanduser("/scrapp2/mtschmitz/")
adata=sc.read(os.path.join(basepath,'macaqueDevBrain-concat_files-MICROGLIA.h5ad'))
sc.settings.figdir="/scrapp2/mtschmitz/Basic_600_4_4_Microglia/"

print([x in adata.var.index for x in ['UTY','SRY','KDM5D','DDX3Y','ZFY']])
pd.DataFrame(adata[:,['UTY','SRY','KDM5D','DDX3Y']].X.todense(),columns=['UTY','SRY','KDM5D','DDX3Y'],index=adata.obs.index)
adata.obs.join(pd.DataFrame(adata[:,['UTY','SRY','KDM5D','DDX3Y']].X.todense(),columns=['UTY','SRY','KDM5D','DDX3Y'],index=adata.obs.index))

adata.obs=adata.obs.join(pd.DataFrame(adata[:,['UTY','SRY','KDM5D','DDX3Y','ZFY']].X.todense(),columns=['UTY','SRY','KDM5D','DDX3Y','ZFY'],index=adata.obs.index))
sc.pl.violin(adata,groupby='timepoint',keys=['UTY','SRY','KDM5D','DDX3Y','ZFY'],save='_seq_genes')

sc.pp.filter_genes(adata,min_counts=1)
adata=sc_utils.std_norm_transform(adata,n_top_genes=8000,tsne=True,save=True)

if 'X_tsne' in adata.obsm.keys():
    sc.pl.tsne(adata,color='batch',save='_batchcolor_MGonly')
    sc.pl.tsne(adata,color='region',save='_regioncolor_MGonly')
    sc.pl.tsne(adata,color='timepoint',save='_timepointcolor_MGonly')
    sc.pl.tsne(adata,color='louvain',save='_louvaincolor_MGonly')

sc.pl.umap(adata,color='batch',save='_batchcolor_MGonly')
sc.pl.umap(adata,color='region',save='_regioncolor_MGonly')
sc.pl.umap(adata,color='timepoint',save='_timepointcolor_MGonly')

adata=sc_utils.cell_cycle_score(adata)
adata=sc_utils.log_reg_diff_exp(adata)
os.rename(os.path.join(sc.settings.figdir,"LouvainLogRegMarkers.csv"),os.path.join(sc.settings.figdir,"PreregressoutLouvainLogRegMarkers.csv"))
os.rename(os.path.join(sc.settings.figdir,"tsne_cc.pdf"),os.path.join(sc.settings.figdir,"tsne_cc_preregress.pdf"))

#adata=sc_utils.marker_analysis(adata,markerpath="~/markers/Markers.txt")

sc.pl.violin(adata,groupby='batch',keys='n_counts',save='ncount')

sc.tl.rank_genes_groups(adata, 'batch', method='logreg')
result = adata.uns['rank_genes_groups']
groups = result['names'].dtype.names
df=pd.DataFrame({group + '_' + key[:1]: result[key][group] for group in groups for key in ['names', 'scores']})
df.to_csv(os.path.join(sc.settings.figdir,"BatchLogRegMarkers.csv"))

sc.pp.regress_out(adata,keys=['batch'])
#sc.pp.scale(adata)
sc.tl.pca(adata)
sc.pp.neighbors(adata, n_neighbors=100)
sc.tl.louvain(adata)
sc.tl.tsne(adata)
sc.pl.tsne(adata,color='batch',save='_batchcolor_MGonly_batchregressed')
sc.pl.tsne(adata,color='louvain',save='_louvaincolor_MGonly_batchregressed')
sc.pl.tsne(adata,color='region',save='_regioncolor_MGonly_batchregressed')
sc.pl.tsne(adata,color='timepoint',save='_timepointcolor_MGonly_batchregressed')
#adata=sc_utils.marker_analysis(adata,markerpath="~/markers/Markers.txt")
adata=sc_utils.cell_cycle_score(adata)
os.rename(os.path.join(sc.settings.figdir,"tsne_cc.pdf"),os.path.join(sc.settings.figdir,"tsne_cc_postregress.pdf"))
adata=sc_utils.log_reg_diff_exp(adata)
os.rename(os.path.join(sc.settings.figdir,"LouvainLogRegMarkers.csv"),os.path.join(sc.settings.figdir,"PostRegressoutLouvainLogRegMarkers.csv"))
sc.tl.phate(adata)
sc.pl.phate(adata)
sc.tl.paga(adata, groups='louvain')
sc.pl.paga(adata, threshold=0.03)
