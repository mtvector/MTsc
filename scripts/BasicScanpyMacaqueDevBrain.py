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

for f in fileList:
    print(f)
    resultList = []
    directory=os.path.join(basepath,'BasicAnalysis',f)
    if not os.path.exists(directory):
        os.makedirs(directory)
    sc.settings.figdir=directory
    sc_filepath=os.path.join(basepath,'Exonic',f,'outs','raw_gene_bc_matrices_h5.h5')
    print(sc_filepath)
    adata=sc_utils.import_file(filename=sc_filepath,refname=refname,save=False)
    print("1")
    adata=sc_utils.quantify_ribo(adata)
    print("2")
    adata=sc_utils.quantify_mito(adata)
    print("3")
    #adata=sc_utils.knee_filter_cells(adata)
    print("4")
    adata=sc_utils.filter_custom(adata,n_counts=400,percent_mito=.6,percent_ribo=.6)
    print(adata)
    adata=sc_utils.sc_lda(adata,n_components=12)
    adata=sc_utils.std_norm_transform(adata,n_top_genes=10000)
    print("5")
    if 'X_tsne' in adata.obsm.keys():
        sc.pl.tsne(adata, color=[x for x in adata.obs.keys() if "lda" in x],save="_lda")
    sc.pl.umap(adata, color=[x for x in adata.obs.keys() if "lda" in x],save="_lda")
    adata=sc_utils.cell_cycle_score(adata)
    print("6")
    adata=sc_utils.log_reg_diff_exp(adata)
    print("7")
    adata=sc_utils.marker_analysis(adata,markerpath="~/markers/Markers.txt",save=True)
    logg.info(os.path.join(basepath,'Exonic',f,'/outs/raw_gene_bc_matrices_h5.h5'))
