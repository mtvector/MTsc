#Script to test cell calling
#Runs a bunch of algos with different quantiles
import os
import scanpy
import scanpy.api as sc
import pandas as pd
import numpy as np
import scipy
from scipy import stats
import re
import loompy
import sklearn
import umi_tools
import matplotlib
import matplotlib.pyplot as plt
matplotlib.use('Agg')
from umi_tools import umi_methods
from collections import Counter
import logging as logg
import importlib.util
spec = importlib.util.spec_from_file_location("ScanpyUtilsMT", os.path.expanduser("~/ScanpyUtilsMT.py"))
sc_utils = importlib.util.module_from_spec(spec)
spec.loader.exec_module(sc_utils)

logg.basicConfig(level=logg.INFO)

refname="refdata-celranger-mmul8-toplevel"

Qlist=[True,False]

sc.settings.autosave=True
sc.settings.autoshow=False

basepath=os.path.expanduser("~/tmp/Macaque/")
fileList = ["E40_MGE_Out","E40_V1_Out","E100temporal_Out_180326_A00269_0056_AH7Y2KDMXX-E100temporal_S16","E100hippo_Out",
           "E40_DRG_Out","E40_Pre-optic_Out","orangutanorganoid_Mmul8_Out"]

fileList=os.listdir(os.path.join(basepath,'Exonic'))
fileList=[os.path.join(basepath,'Exonic',f,'outs','raw_gene_bc_matrices_h5.h5') for f in fileList]


resultList = []
directory=os.path.join(basepath,'BasicAnalysisAggregate_600_4_4')
if not os.path.exists(directory):
    os.makedirs(directory)
sc.settings.figdir=directory

if not os.path.exists(os.path.join(directory,'macaqueDevBrain-concat_files-do_pca-do_tsne-neighbors_and_umap-louvain.h5ad')):
    adata=sc_utils.concat_files(fileList=fileList,refname=refname,groupname="macaqueDevBrain",n_counts=600,percent_mito=.4,percent_ribo=.4,filter_genes=10,save=True)
    print("1")
    sc.pp.filter_genes(adata,min_counts=1)
    print(adata)
    adata=sc_utils.sc_lda(adata,n_components=12)
    print("2")
    adata=sc_utils.std_norm_transform(adata,n_top_genes=6000,tsne=True,save=True)
    print("5")
else:
    adata=sc.read(os.path.join(directory,'macaqueDevBrain-concat_files-do_pca-do_tsne-neighbors_and_umap-louvain.h5ad'))
    print(adata)

adata.obs.region=[x.lower() for x in adata.obs.region]

if 'X_tsne' in adata.obsm.keys():
    sc.pl.tsne(adata,color='batch',save='batchcolor')
    sc.pl.tsne(adata,color='region',save='regioncolor')
    sc.pl.tsne(adata,color='timepoint',save='timepointcolor')

sc.pl.umap(adata,color='batch',save='batchcolor')
sc.pl.umap(adata,color='region',save='regioncolor')
sc.pl.umap(adata,color='timepoint',save='timepointcolor')


if 'X_tsne' in adata.obsm.keys():
    sc.pl.tsne(adata, color=[x for x in adata.obs.keys() if "lda" in x],save="_lda")
sc.pl.umap(adata, color=[x for x in adata.obs.keys() if "lda" in x],save="_lda")

adata=sc_utils.cell_cycle_score(adata)
print("6")
adata=sc_utils.log_reg_diff_exp(adata)
print("7")
adata=sc_utils.marker_analysis(adata,markerpath="~/markers/Markers.txt")
