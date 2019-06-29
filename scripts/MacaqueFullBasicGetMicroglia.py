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

the_filename=os.path.expanduser('~/subsetCells80k.txt')
if not os.path.exists(the_filename):
    subset = random.sample(adata.obs.index.tolist(),6000)
    with open(the_filename, 'w') as f:
        for s in subset:
            f.write(s + '\n')
else:
    with open(the_filename, 'r') as f:
        subset = [line.rstrip('\n') for line in f]


#Create a list of the raw matrices I'm going to be looking at
filepath=os.path.expanduser("/ye/yelabstore2/mtschmitz/seq/AlignedMacaqueBrainSeq/")
fileNameList=os.listdir(os.path.join(filepath,'Exonic'))
fileList=[os.path.join(filepath,'Exonic',f,'outs','raw_gene_bc_matrices_h5.h5') for f in fileNameList]

#Import and concatenate whole list of files into a single unified object with correct metadata (carries out some minor filtering too, see function)
basepath=os.path.expanduser("/scrapp2/mtschmitz/")
adataFull=sc.read(os.path.join(basepath,'macaqueDevBrain-concat_files.h5ad'))
#adataFull=adataFull[subset,:]
#adataFull.write(os.path.join(basepath,'macaqueDevBrain-concat_files_SUBSET.h5ad'))

#Set the directory to the place we want to save outputs
directory=os.path.join(basepath,'Basic_600_4_4')
adata=sc.read(os.path.join(directory,'macaqueDevBrain-concat_files-do_pca-do_tsne-neighbors_and_umap-louvain.h5ad'))
#adata=adata[adataFull.obs.index.tolist(),:]
adata=sc_utils.marker_analysis(adata,markerpath='~/markers/Markers.txt',save=True)
#adata.write(os.path.join(directory,'macaqueDevBrain-concat_files-do_pca-do_tsne-neighbors_and_umap-louvain_SUBSET.h5ad'))
directory=os.path.join(basepath,'Basic_600_4_4_Microglia')
sc.settings.figdir=directory
#sc.pl.umap(adata,color=["louvain","percent_mito",'percent_ribo',"n_counts"], save="_louvain")
#sc.pl.umap(adata,color='batch',save='_batchcolor')
#sc.pl.umap(adata,color='region',save='_regioncolor')
#sc.pl.umap(adata,color='timepoint',save='_timepointcolor')

adataFull.obs.region=[x.lower() for x in adataFull.obs.region]
#Simple gene filtering (speeds up norm and transform a little)
sc.pp.filter_genes(adataFull,min_cells=15)
#save the raw data
adataFull.raw=adata
sc.settings.figdir=directory

if not os.path.exists(directory):
    os.makedirs(directory)
grouped=adata.obs.groupby(['louvain']).mean()
print(grouped)
grouped=grouped.loc[:,adata.uns['marker_groups']]
print(grouped)
print(np.array(adata.obs['louvain']))
print(np.array(grouped.idxmax(axis=0)))
print(np.array(grouped.idxmax(axis=0)=='Microglia'))
print(np.array(grouped.idxmax(axis=1)))
print(np.array(grouped.idxmax(axis=1)=='Microglia'))
print(grouped.index[np.array(grouped.idxmax(axis=1)=='Microglia')].tolist())
sys.stdout.flush()
cellIndexList=grouped.index.values[np.array(grouped.idxmax(axis=1)=='Microglia')].tolist()
cell_clusters=np.array([x in cellIndexList for x in adata.obs['louvain']])
print(adata.obs.index.values[cell_clusters])
sys.stdout.flush()
#adataFull=adataFull[adata.obs.index.values[cell_clusters],:]
#adata=adata[adata.obs.index.values[cell_clusters],:]
adataFull._inplace_subset_obs(adata.obs.index.values[cell_clusters])
print(adataFull)
adataFull.write(os.path.join(basepath,'macaqueDevBrain-concat_files-MICROGLIA.h5ad'))
adata._inplace_subset_obs(adata.obs.index.values[cell_clusters])
print('written')
sys.stdout.flush()
sc.pl.umap(adata,color=["louvain","percent_mito",'percent_ribo',"n_counts"], save="_louvain_micro")
sc.pl.umap(adata,color='batch',save='_batchcolor_micro')
sc.pl.umap(adata,color='region',save='_regioncolor_micro')
sc.pl.umap(adata,color='timepoint',save='_timepointcolor_micro')
print('written')
sys.stdout.flush()
sc.pl.umap(adata,color=["louvain","percent_mito",'percent_ribo',"n_counts"], save="_louvain_micro")
sc.pl.umap(adata,color='batch',save='_batchcolor_micro')
sc.pl.umap(adata,color='region',save='_regioncolor_micro')
sc.pl.umap(adata,color='timepoint',save='_timepointcolor_micro')
