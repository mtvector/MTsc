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

filepath=os.path.expanduser("/ye/yelabstore2/mtschmitz/seq/AlignedMacaqueBrainSeq/")
fileNameList=os.listdir(os.path.join(filepath,'Exonic'))
fileList=[os.path.join(filepath,'Exonic',f,'outs','raw_gene_bc_matrices_h5.h5') for f in fileNameList]

#Import and concatenate whole list of files into a single unified object with correct metadata (carries out some minor filtering too, see function)
basepath=os.path.expanduser("/scrapp2/mtschmitz/")

resultList = []
directory=os.path.join(basepath,'BasicAnalysisAggregateMNN_500_5_5')
if not os.path.exists(directory):
    os.makedirs(directory)
sc.settings.figdir=directory

adata=sc_utils.concat_files_mnn(fileList=fileList[0:3],refname=refname,n_top_genes=6000,groupname="macaqueDevBrain",n_counts=800,percent_mito=.4,percent_ribo=.4,filter_genes=15,save=True)
adata=sc_utils.std_norm_transform(adata,n_top_genes=3000,log=False)
sc.pl.tsne(adata,color='batch',save='smallbatchcolor')

adata,mnn_list,angle_list=sc_utils.concat_files_mnn(fileList=fileList,refname=refname,n_top_genes=3000,groupname="macaqueDevBrain",n_counts=800,percent_mito=.4,percent_ribo=.4,filter_genes=15,save=True)
print("1")
print(adata)
print(mnn_list)
print(angle_list)
adata=sc_utils.std_norm_transform(adata,n_top_genes=3000,log=False)
sc.pl.tsne(adata,color='batch',save='batchcolor')
print("5")
adata=sc_utils.cell_cycle_score(adata)
print("6")
adata=sc_utils.marker_analysis(adata,markerpath="~/markers/Markers.txt")
print("7")
adata=sc_utils.log_reg_diff_exp(adata,save=True)
