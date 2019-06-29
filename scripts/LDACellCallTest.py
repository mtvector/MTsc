import os
import scanpy
import scanpy.api as sc
import pandas as pd
import numpy as np
import sklearn
from sklearn import decomposition
import matplotlib
import matplotlib.pyplot as plt
matplotlib.use('Agg')
plt.switch_backend('agg')
import logging as logg
import seaborn


logg.basicConfig(level=logg.INFO)

def adataToTrainState(filename,refname,num_genes=8000):
    adata = sc.read_10x_h5(filename,refname)
    adata.name=filename.split(os.sep)[-3]
    sc.tl.addCleanObsNames(adata)
    sc.pp.filter_cells(adata, min_genes=10)
    mito_genes = [name for name in adata.var_names if name in ['ND1','ND2','ND4L','ND4','ND5','ND6','ATP6','ATP8','CYTB','COX1','COX2','COX3']]
    #mito_genes = [name for name in adata.var_names if name.startswith(('MTND','MTCO','MTATP','MTCYB','MTRNR','MTTR',))]
    ribo_genes = [name for name in adata.var_names if name.startswith('RPS') or name.startswith('RPL') ]
    adata.obs['percent_mito'] = np.sum(
        adata[:, mito_genes].X, axis=1).A1 / np.sum(adata.X, axis=1).A1
    adata.obs['percent_ribo'] = np.sum(
        adata[:, ribo_genes].X, axis=1).A1 / np.sum(adata.X, axis=1).A1
    # add the total counts per cell as observations-annotation to adata
    adata.obs['n_counts'] = adata.X.sum(axis=1).A1
    scanpy.api.pp.filter_genes_dispersion(adata,n_top_genes=min(np.sum(np.sum(adata.X, axis=0)>0),num_genes))
    return(adata)

def moving_average(data,window_width):
    cumsum_vec = np.cumsum(np.insert(data, 0, 0))
    ma_vec = (cumsum_vec[window_width:] - cumsum_vec[:-window_width]) / window_width
    return(ma_vec)

refname="refdata-celranger-mmul8-toplevel"

sc.settings.autosave=True
sc.settings.autoshow=False

basepath=os.path.expanduser("~/tmp/Macaque/")
fileList=os.listdir(os.path.join(basepath,'Exonic'))

for f in fileList:
    resultList = []
    directory=os.path.join(basepath,'CellCallLDA',f)
    adata=adataToTrainState(os.path.join(basepath,'Exonic',f,'outs/raw_gene_bc_matrices_h5.h5'),refname)
    sc.settings.figdir=directory
    if not os.path.exists(directory):
        os.makedirs(directory)
    print(os.path.join(basepath,'Exonic',f,'/outs/raw_gene_bc_matrices_h5.h5'))
    name=adata.name
    ldaM=sklearn.decomposition.LatentDirichletAllocation(n_components=10,learning_method="online",verbose=2)
    ldaM.fit(adata.X)
    doc_topic=ldaM.transform(adata[adata.obs.n_counts.argsort(),:].X)
    hm=seaborn.heatmap(doc_topic)
    hm.savefig(os.path.join(sc.settings.figdir,name+"_TopicHeatmap.png"))
    plt.plot(range(doc_topic.shape[0]),np.std( doc_topic,axis=1))
    plt.title("SD of topic allocation in cells")
    plt.xlabel("RankedCells")
    plt.ylabel("SD")
    plt.savefig(os.path.join(sc.settings.figdir,name+"_TopicCellSD"))

    convolvedSD=moving_average(np.std( doc_topic,axis=1),500)
    plt.plot(range(len(convolvedSD)),convolvedSD)
    plt.title("Moving Average Sparsity (SD)")
    plt.xlabel("RankedCells")
    plt.ylabel("SD")
    plt.savefig(os.path.join(sc.settings.figdir,name+"_TopicCellMASD.png"))
    plt.plot(range(adata.obs.shape[0]),np.sort(np.log10(adata.obs.n_counts))[::-1])
    plt.title("Log Counts per Cell")
    plt.xlabel("RankedCells")
    plt.title("Log Reads")
    plt.savefig(os.path.join(sc.settings.figdir,name+"_LogCounts.png"))

    for i in range(doc_topic.shape[1]):
        convolved=moving_average(doc_topic[:,i],500)
        plt.plot(range(len(convolved)),convolved)
    plt.xlabel("RankedCells")
    plt.ylabel("Moving Avg Topic Probability")
    plt.savefig(os.path.join(sc.settings.figdir,name+"_TopicMovingAvg.png"))
