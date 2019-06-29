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
import ruptures as rpt
import pickle

logg.basicConfig(level=logg.INFO)

def adataToTrainStateAggregate(fileList,refname,num_genes=3000):
    adatas=[]
    batch_categories=[]
    for filename in fileList:
        a = sc.read_10x_h5(filename,refname)
        sc.pp.filter_cells(a, min_genes=30)
        a.name=filename.split(os.sep)[-3]
        a.obs['dataset']=a.name
        batch_categories.append(a.name)
        adatas.append(a)
    adata = sc.AnnData.concatenate(*adatas, batch_categories=batch_categories)
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
    return(adata[adata.obs.n_counts.argsort()[::-1],:])

def moving_average(data,window_width):
    cumsum_vec = np.cumsum(np.insert(data, 0, 0))
    ma_vec = (cumsum_vec[window_width:] - cumsum_vec[:-window_width]) / window_width
    return(ma_vec)

def print_top_words(model, feature_names, n_top_words):
    for topic_idx, topic in enumerate(model.components_):
        message = "Topic #%d: " % topic_idx
        message += " ".join([feature_names[i]
                             for i in topic.argsort()[:-n_top_words - 1:-1]])
        print(message)
    print()

def table_top_words(model, feature_names, n_top_words):
    df=[]
    for topic_idx, topic in enumerate(model.components_):
        print(topic_idx)
        df.append([feature_names[i]
                            for i in topic.argsort()[:-n_top_words - 1:-1]])
    print(df)
    return(pd.DataFrame(df))

refname="refdata-celranger-mmul8-toplevel"

sc.settings.autosave=True
sc.settings.autoshow=False

basepath=os.path.expanduser("~/tmp/Macaque/")
nameList=os.listdir(os.path.join(basepath,'Exonic'))
fileList=[os.path.join(basepath,'Exonic',f,'outs/raw_gene_bc_matrices_h5.h5') for f in nameList]

f="MacaqueDevBrain"
resultList = []
directory=os.path.join(basepath,'AggregateCellCallLDA',f)
adata=adataToTrainStateAggregate(fileList,refname)
print(adata)
sc.settings.figdir=directory
if not os.path.exists(directory):
    os.makedirs(directory)
name="MacaqueBrainDev"
ldaM=sklearn.decomposition.LatentDirichletAllocation(n_components=10,learning_method="online",verbose=2,n_jobs=-1)
ldaM.fit(adata.X)
print("Model done")
pickle.dump(ldaM, open(os.path.join(directory,"LDAmodel"), 'wb'))
doc_topic=ldaM.transform(adata.X)
hm=seaborn.heatmap(doc_topic).get_figure()
hm.savefig(os.path.join(sc.settings.figdir,name+"_TopicHeatmap.png"))
plt.clf()
plt.plot(range(doc_topic.shape[0]),np.std( doc_topic,axis=1))
plt.title("SD of topic allocation in cells")
plt.xlabel("RankedCells")
plt.ylabel("SD")
plt.savefig(os.path.join(sc.settings.figdir,name+"_TopicCellSD"))
plt.clf()
convolvedSD=moving_average(np.std( doc_topic,axis=1),300)
plt.plot(range(len(convolvedSD)),convolvedSD)
plt.title("Moving Average Sparsity (SD)")
plt.xlabel("RankedCells")
plt.ylabel("SD")
plt.savefig(os.path.join(sc.settings.figdir,name+"_TopicCellMASD.png"))
plt.clf()
plt.plot(range(adata.obs.shape[0]),np.log10(adata.obs.n_counts))
plt.title("Log Counts per Cell")
plt.xlabel("RankedCells")
plt.title("Log Reads")
plt.savefig(os.path.join(sc.settings.figdir,name+"_LogCounts.png"))
plt.clf()
for i in range(doc_topic.shape[1]):
    convolved=moving_average(doc_topic[:,i],300)
    plt.plot(range(len(convolved)),convolved)
plt.xlabel("RankedCells")
plt.ylabel("Moving Avg Topic Probability")
plt.savefig(os.path.join(sc.settings.figdir,name+"_TopicMovingAvg.png"))
plt.clf()
convolvedSD=moving_average(adata.obs['percent_ribo'].tolist(),300)
plt.plot(range(len(convolvedSD)),convolvedSD)
plt.title("Moving Average Percent Ribo")
plt.savefig(os.path.join(sc.settings.figdir,name+"_RiboCounts.png"))
plt.clf()
convolvedSD=moving_average(adata.obs['percent_mito'].tolist(),300)
plt.plot(range(len(convolvedSD)),convolvedSD)
plt.title("Moving Average Percent Mito")
plt.savefig(os.path.join(sc.settings.figdir,name+"_MitoCounts.png"))
plt.clf()
signal = np.column_stack((np.std( doc_topic,axis=1),
                      adata.obs['percent_ribo'].tolist(),
                      adata.obs['percent_mito'].tolist(),
                      np.array(list(adata.obs.n_counts))))

algo = rpt.Window(width=2500,model="l1").fit(signal)
result = algo.predict(pen=50)
costs=[]
for i in range(len(result)-1):
    costs.append(algo.cost.sum_of_costs([result[i],result[len(result)-1]]))
rpt.display(signal=signal,true_chg_pts=result,computed_chg_pts= result)
plt.title(str(np.argmin(costs))+' <- Best cPoint')
plt.savefig(os.path.join(sc.settings.figdir,name+"_Changepoints.png"))
plt.clf()
print_top_words(ldaM,adata.var.index,15)
table_top_words(ldaM,adata.var.index,25).to_csv(os.path.join(sc.settings.figdir,name+"_TopicMarkers.txt"))
