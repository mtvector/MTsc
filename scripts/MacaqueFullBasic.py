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
import importlib.util
spec = importlib.util.spec_from_file_location("ScanpyUtilsMT", os.path.expanduser("~/ScanpyUtilsMT.py"))
sc_utils = importlib.util.module_from_spec(spec)
spec.loader.exec_module(sc_utils)

logg.basicConfig(level=logg.INFO)
sc.settings.verbosity=5

refname="refdata-celranger-mmul8-toplevel"

sc.settings.autosave=True
sc.settings.autoshow=False

filepath=os.path.expanduser("/ye/yelabstore2/mtschmitz/seq/AlignedMacaqueBrainSeq/")
fileNameList=os.listdir(os.path.join(filepath,'Exonic'))
fileList=[os.path.join(filepath,'Exonic',f,'outs','raw_gene_bc_matrices_h5.h5') for f in fileNameList]


basepath=os.path.expanduser("/scrapp2/mtschmitz/")
sc.settings.figdir=os.path.expanduser('~/')
if not os.path.exists(os.path.join(basepath,'macaqueDevBrain-concat_files.h5ad')):
    adata=sc_utils.concat_files(fileList=fileList,refname=refname,groupname="macaqueDevBrain",n_counts=600,percent_mito=.4,percent_ribo=.4,filter_genes=10,save=True)
    import shutil
    shutil.copyfile(os.path.expanduser('~/macaqueDevBrain-concat_files.h5ad'),os.path.join(basepath,'macaqueDevBrain-concat_files.h5ad'))
else:
    adata=sc.read(os.path.join(basepath,'macaqueDevBrain-concat_files.h5ad'))


directory=os.path.join(basepath,'Basic_600_4_4')
if not os.path.exists(directory):
    os.makedirs(directory)
sc.settings.figdir=directory
adata=sc.read(os.path.join(basepath,'macaqueDevBrain-concat_files.h5ad'))
adata.obs.region=[x.lower() for x in adata.obs.region]
sc.pp.filter_genes(adata,min_cells=15)
adata.raw=adata
adata=sc_utils.std_norm_transform(adata,n_top_genes=6000,tsne=True,save=True)
if 'X_tsne' in adata.obsm.keys():
    sc.pl.tsne(adata,color=["louvain","percent_mito",'percent_ribo',"n_counts"], save="_louvain")
    sc.pl.tsne(adata,color='batch',save='_batchcolor')
    sc.pl.tsne(adata,color='region',save='_regioncolor')
    sc.pl.tsne(adata,color='timepoint',save='timepointcolor')
sc.pl.umap(adata,color=["louvain","percent_mito",'percent_ribo',"n_counts"], save="_louvain")
sc.pl.umap(adata,color='batch',save='_batchcolor')
sc.pl.umap(adata,color='region',save='_regioncolor')
sc.pl.umap(adata,color='timepoint',save='_timepointcolor')
a=a=adata.obs.groupby(['louvain','region']).size().unstack(fill_value=0)
hm=seaborn.heatmap((a.T/a.sum(axis=1)).T).get_figure()
hm.savefig(os.path.join(sc.settings.figdir,"RegionCluster.png"))
plt.clf()
pd.DataFrame(adata.obs.groupby(['louvain','region']).size().unstack(fill_value=0)).to_csv(os.path.join(sc.settings.figdir,"RegionCluster.csv"))
pd.DataFrame(adata.obs.groupby(['louvain','batch']).size().unstack(fill_value=0)).to_csv(os.path.join(sc.settings.figdir,"BatchCluster.csv"))
adata=sc_utils.marker_analysis(adata,markerpath='~/markers/Markers.txt',save=False)
adata=sc_utils.log_reg_diff_exp(adata,save=True)
