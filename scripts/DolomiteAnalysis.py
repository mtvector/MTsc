#Script to test cell calling
#Runs a bunch of algos with different quantiles
import os
import scanpy
import scanpy.api as sc
import pandas as pd
import numpy as np
import re
import loompy
import matplotlib
import matplotlib.pyplot as plt
matplotlib.use('Agg')
from collections import Counter
import logging as logg
import importlib.util
spec = importlib.util.spec_from_file_location("ScanpyUtilsMT", os.path.expanduser("~/utils/ScanpyUtilsMT.py"))
sc_utils = importlib.util.module_from_spec(spec)
spec.loader.exec_module(sc_utils)
logg.basicConfig(level=logg.INFO)


sc.settings.autosave=True
sc.settings.autoshow=False

basepath=os.path.expanduser("/home/mschmitz1/matthew/dolo/")

fileList=os.listdir(os.path.join(basepath,'pollena'))
fileList=[s for s in fileList if '_1' in s or 'merged' in s]

for f in fileList:
    print(f)
    resultList = []
    directory=os.path.join(basepath,'BasicAnalysis',f)
    if not os.path.exists(directory):
        os.makedirs(directory)
    sc.settings.figdir=directory
    sc_filepath=os.path.join(basepath,'pollena',f)
    print(sc_filepath)
    adata=sc_utils.import_dropseq(filename=sc_filepath,filter_genes=500,save=False)
    print("1")
    adata=sc_utils.quantify_ribo(adata)
    print("2")
    adata=sc_utils.quantify_mito(adata)
    print("3")
    #adata=sc_utils.knee_filter_cells(adata)
    print("4")
    adata=sc_utils.filter_custom(adata,n_counts=600,percent_mito=.6,percent_ribo=.6)
    adata=sc_utils.std_norm_transform(adata,tsne=True,n_top_genes=10000)
    sc.pl.violin(adata, ['n_genes', 'n_counts', 'percent_mito','percent_ribo'],
             jitter=0.4, multi_panel=True,save='stats')
    sc.pl.scatter(adata, x='n_counts', y='percent_mito',save='percentMito')
    sc.pl.scatter(adata, x='n_counts', y='n_genes',save='numgenes')
    print("5")
    print(adata.var)
    print(adata.obs)
    print(adata.X)
    adata=sc_utils.cell_cycle_score(adata)
    print("6")
    if len(set(adata.obs.leiden))>1:
        adata=sc_utils.log_reg_diff_exp(adata)
    print("7")
    adata=sc_utils.marker_analysis(adata,variables=['leiden'],markerpath="~/markers/Markers.txt",save=False)
    adata.obs.to_csv(os.path.join(sc.settings.figdir,'ObsTable.csv'))
