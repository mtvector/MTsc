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
print("start")
logg.basicConfig(level=logg.INFO)
refname="refdata-celranger-mmul8-toplevel"

sc.settings.autosave=True
sc.settings.autoshow=False
sc.settings.verbosity=5

filepath=os.path.expanduser("/ye/yelabstore2/mtschmitz/seq/AlignedMacaqueBrainSeq/")
fileList=os.listdir(os.path.join(filepath,'Exonic'))
fileList=[os.path.join(filepath,'Exonic',f,'outs','raw_gene_bc_matrices_h5.h5') for f in fileList]

print('ReadingAdata')
basepath=os.path.expanduser("/scrapp2/mtschmitz/HDPh5ad/")
sc.settings.figdir=os.path.expanduser('~/')
if not os.path.exists(os.path.join(basepath,'macaqueDevBrain-concat_files.h5ad')):
    adata=sc_utils.concat_files(fileList=fileList,refname=refname,groupname="macaqueDevBrain",n_counts=600,percent_mito=.4,percent_ribo=.4,filter_genes=10,save=True)
    import shutil
    shutil.copyfile(os.path.expanduser('~/macaqueDevBrain-concat_files.h5ad'),os.path.join(basepath,'macaqueDevBrain-concat_files.h5ad'))

#Read in subset of cells to use
the_filename=os.path.expanduser('~/subsetCells80k.txt')
if not os.path.exists(the_filename):
    subset = random.sample(adata.obs.index.tolist(),80000)
    with open(the_filename, 'w') as f:
        for s in subset:
            f.write(s + '\n')
else:
    with open(the_filename, 'r') as f:
        subset = [line.rstrip('\n') for line in f]

print('StartLoop')
sys.stdout.flush()
'''
#Loop to check
defaults=[1,1,0]
v=[0.01,1,100]
test_condit=np.empty([len(v)*len(defaults),len(v)])
for i in range(len(v)):
    for j in range(len(defaults)):
        test_condit[(i*len(v)):((i+1)*len(v)),j]= v if i==j else  np.repeat(repeats=len(v),a=v[defaults[j]])
'''
test_condit=[[1,.1,10],[10,.1,1],[10,.1,10],[1,1,1],[1,1,.01],[.01,.01,.01],[100,100,.01],[1,.1,1],[1,1,10],[1,1,.001],[1,1,1],[.1,1,.01],[1,.1,.01],[10,1,.01],[1,10,.01]]

for alpha,gamma,eta in test_condit:
    directory=os.path.join(basepath,'HDPBasicAnalysisAggregateDownsample_600_4_4_'+str(alpha)+'alpha_'+str(gamma)+'gamma_'+str(eta)+'eta')
    if not os.path.exists(directory):
        os.makedirs(directory)
    sc.settings.figdir=directory
    if not os.path.exists(os.path.join(directory,'macaqueDevBrain-concat_files-sc_hdp.h5ad')):
        adata=sc.read(os.path.join(basepath,'macaqueDevBrain-concat_files.h5ad'))
        adata=adata[subset,:]
        adata.obs.region=[x.lower() for x in adata.obs.region]
        sc.pp.filter_genes(adata,min_cells=15)
        adata=sc_utils.sc_hdp(adata,alpha=alpha,eta=eta,gamma=gamma,eps=1e-5,save=True)
    else:
        adata=sc.read(os.path.join(directory,'macaqueDevBrain-concat_files-sc_hdp.h5ad'))
    adata.raw=adata
    sc_utils.dirichlet_marker_analysis(adata,markerpath='~/markers/Markers.txt')
    print(adata)
    logg.info('KNNING!')
    print(adata)
    sys.stdout.flush()
    neighbors=sc.Neighbors(anndata.AnnData(adata.obsm['cell_topic']))
    neighbors.compute_neighbors(n_neighbors=100,use_rep='X')
    print('KNNed')
    sys.stdout.flush()
    adata.uns['neighbors'] = {}
    adata.uns['neighbors']['distances'] = neighbors.distances
    adata.uns['neighbors']['connectivities'] = neighbors.connectivities
    print('louvaining')
    sys.stdout.flush()
    sc.tl.louvain(adata)
    pd.DataFrame(adata.obs.groupby(['louvain','region']).size().unstack(fill_value=0)).to_csv(os.path.join(sc.settings.figdir,"RegionCluster.csv"))
    print('umapping')
    sc.tl.umap(adata)
    sc.pl.umap(adata,color=["louvain"], save="_HDPKNN_louvain")
    sc.pl.umap(adata,color=["louvain","percent_mito",'percent_ribo',"n_counts"], save="_HDPKNN_stats")
    sc.pl.umap(adata,color='batch',save='_HDPKNN_batchcolor')
    sc.pl.umap(adata,color='region',save='_HDPKNN_regioncolor')
    sc.pl.umap(adata,color='timepoint',save='_HDPKNN_timepointcolor')
    sc.pp.normalize_per_cell(adata, key_n_counts='n_counts')
    sc.pp.log1p(adata)
    sc.pp.scale(adata)
    print(adata)
    sc_utils.marker_analysis(adata,markerpath='~/markers/Markers.txt',save=False)
    a=a=adata.obs.groupby(['louvain','region']).size().unstack(fill_value=0)
    hm=seaborn.heatmap((a.T/a.sum(axis=1)).T).get_figure()
    hm.savefig(os.path.join(sc.settings.figdir,"RegionCluster.png"))
    plt.clf()
    pd.DataFrame(adata.obs.groupby(['louvain','region']).size().unstack(fill_value=0)).to_csv(os.path.join(sc.settings.figdir,"RegionCluster.csv"))
    pd.DataFrame(adata.obs.groupby(['louvain','batch']).size().unstack(fill_value=0)).to_csv(os.path.join(sc.settings.figdir,"BatchCluster.csv"))
    a=a=adata.obs.groupby(['louvain','batch']).size().unstack(fill_value=0)
    hm=seaborn.heatmap((a.T/a.sum(axis=1)).T).get_figure()
    hm.savefig(os.path.join(sc.settings.figdir,"BatchCluster.png"))
    plt.clf()
    #sc_utils.log_reg_diff_exp(adata)
