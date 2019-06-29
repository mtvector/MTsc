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
#adata=sc.read(os.path.join('/scrapp2/mtschmitz/HDPh5ad/HDPBasicAnalysisAggregateDownsample_600_4_4_0.1alpha_1.0gamma_0.1eta/macaqueDevBrain-concat_files-sc_hdp.h5ad'))

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


directory=os.path.join(basepath,'SplitTimeBasic_600_4_4_E40.0'))
sc.settings.figdir=directory
adata=sc.read(os.path.join(basepath,'macaqueDevBrain-concat_files-do_pca-do_tsne-neighbors_and_umap-louvain-marker_analysis-log_reg_diff_exp.h5ad'))
sc.pl.tsne(adata,color=["louvain",'region','PPP1R1B',"ARHGAP11A",'RAB13','ZEB1','PTPN11','SPDL1',"PAX8","PAX3","GNG8","RERE","APOE"], save="_checkmarkers")
