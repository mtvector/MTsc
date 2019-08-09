import scanpy as sc
import scvelo as scv
import anndata
import os
import re
import pandas as pd
import numpy as np
import scipy
import matplotlib.pyplot as plt
import seaborn as sns
#scv.settings.set_figure_params('scvelo')
sc.settings.figdir=str(os.path.expanduser('~/figs/'+"MultiseqSCvelo"))
sc.settings.autosave=True
sc.settings.autoshow=False
import scvelo as scv
#scv.settings.figdir=sc.settings.figdir
scv.settings.autosave=True
scv.settings.autoshow=False
if not os.path.exists(sc.settings.figdir):
        os.makedirs(sc.settings.figdir)
        
import importlib.util
spec = importlib.util.spec_from_file_location("ScanpyUtilsMT", os.path.expanduser("~/code/pollye/MTsc/utils/ScanpyUtilsMT.py"))
sc_utils = importlib.util.module_from_spec(spec)
spec.loader.exec_module(sc_utils)
        
headpath=os.path.expanduser("/scrapp2/mtschmitz")

adata=sc.read_h5ad(os.path.join(headpath,'MultiseqRaw.h5ad'))
print(adata)

adata.var_names_make_unique()
adata._inplace_subset_obs([x is not None for x in [re.search('[a-zA-Z]', x) for x in adata.obs['region']]])
adata.obs['region']=[re.sub('_A_|_B_','_',x) for x in adata.obs['region']]
adata._inplace_subset_obs([x is None for x in [re.search('nan|Doublet|Negative', x) for x in adata.obs['region']]])
adata.obs['ActualRegion']=sc_utils.macaque_correct_regions(adata.obs['ActualRegion'])


mito_genes = [name for name in adata.var_names if name in ['ND1','ND2','ND4L','ND4','ND5','ND6','ATP6','ATP8','CYTB','COX1','COX2','COX3'] or 'MT-' in name]
ribo_genes = [name for name in adata.var_names if name.startswith('RPS') or name.startswith('RPL') ]

adata.obs['percent_mito'] = np.sum(
    adata[:, mito_genes].X, axis=1).A1 / np.sum(adata.X, axis=1).A1
adata.obs['percent_ribo'] = np.sum(
    adata[:, ribo_genes].X, axis=1).A1 / np.sum(adata.X, axis=1).A1
# add the total counts per cell as observations-annotation to adata
adata.obs['n_counts'] = adata.X.sum(axis=1).A1
#scanpy.api.pp.filter_genes_dispersion(adata,n_top_genes=np.sum(np.sum(adata.X, axis=0)>0))
sc.pp.filter_genes(adata, min_cells=20,inplace=True)
sc.pp.filter_cells(adata,min_genes=400)
#sc.pl.violin(adata, ['n_genes', 'n_counts', 'percent_mito'],jitter=0.4, multi_panel=True)

sc.pp.normalize_per_cell(adata, counts_per_cell_after=1e4)

sc.pp.log1p(adata)
normalizedadata=adata.copy()
sc.pp.highly_variable_genes(adata, n_top_genes=15000)
sc.pp.scale(adata, max_value=10)
sc.pp.regress_out(adata,'batch')

sc.pp.pca(adata)
sc.pp.neighbors(adata)
sc.tl.umap(adata)
sc.tl.leiden(adata)

sc_utils.cell_cycle_score(adata)
'''
sc.pl.umap(adata, color=['leiden','n_counts','region','phase'])
sc.pl.umap(adata, color=['region','simpleregion','batch','leiden'],save='_MultiseqSummary')
sc.pl.umap(adata, color=['region'],save='_region')
sc.pl.umap(adata, color=['leiden'],save='_leiden')
sc.pl.umap(adata, color=['tp'],save='_tp')
sc.pl.umap(adata, color=['ActualRegion'],save='_ActualRegion')


sc.pl.umap(adata[['Negative' not in x and 'Doublet' not in x for x in adata.obs['region']],:], color=['leiden'],save='_screenedleiden')
sc.pl.umap(adata[['Negative' not in x and 'Doublet' not in x for x in adata.obs['region']],:], color=['region'],save='_screenedsample')
sc.pl.umap(adata[['Negative' not in x and 'Doublet' not in x for x in adata.obs['region']],:], color=['simpleregion'],save='_screenedsimpleregion')
'''
#df=pd.DataFrame([sc_utils.get_cluster_metrics(adata,rands=['batch','ActualRegion','region','tp'])])
#df.to_csv(os.path.join(sc.settings.figdir,"Metrics.csv"))

adata=sc_utils.marker_analysis(adata,variables=['leiden','region'],markerpath=os.path.expanduser('~/markers.txt'),subclass=True)

adata.obs['subclassname']=[re.sub('mkrscore','',x) for x in adata.obs.loc[:,['mkrscore' in x for x in adata.obs.columns]].astype('float').idxmax(axis=1)]
normalizedadata.obs['subclassname']=[re.sub('mkrscore','',x) for x in adata.obs.loc[:,['mkrscore' in x for x in adata.obs.columns]].astype('float').idxmax(axis=1)]
sc.pl.umap(adata,color=['subclassname'],save='_subclassname')

def most_frequent(List): 
    return max(set(List), key = List.count) 
classlist=[]
for c in adata.obs['subclassname']:
    fullbool=[c in x for x in adata.uns['markers']['fullclass']]
    flatclass=[item for sublist in adata.uns['markers'].loc[fullbool,'type'] for item in sublist]
    classlist.append(most_frequent(flatclass))
adata.obs['classname']=classlist
normalizedadata.obs['classname']=classlist
sc.pl.umap(adata,color=['classname'],save='_classname')

def variance(X,axis=1):
    return( (X.power(2)).mean(axis=axis)-(np.power(X.mean(axis=axis),2)) )

key1='region'
key2='classname'
meanmat=np.zeros((len(set(normalizedadata.obs[key1])),len(set(normalizedadata.obs[key2])),normalizedadata.shape[1]))
varmat=np.zeros((len(set(normalizedadata.obs[key1])),len(set(normalizedadata.obs[key2])),normalizedadata.shape[1]))
for ri,r in enumerate(set(normalizedadata.obs[key1])):
    for scni,scn in enumerate(set(adata.obs[key2])):
        print(scn)
        a=normalizedadata[(normalizedadata.obs[key1]==r) & (normalizedadata.obs[key2]==scn),:]
        if a.shape[0]>5:
            meanmat[ri,scni]=variance(a.X,axis=0)/a.X.mean(axis=0)
            varmat[ri,scni]=sc.pp.highly_variable_genes(a,inplace=False)['dispersions_norm']
                       
fig, ax = plt.subplots(figsize=(20,20)) 
sns.heatmap(np.nanmean(varmat,axis=2),yticklabels=set(normalizedadata.obs[key1]),xticklabels=set(normalizedadata.obs[key2]))
plt.savefig(os.path.join(sc.settings.figdir,key1+key2+'DispersionHeatmap.pdf'))
plt.close()

fig, ax = plt.subplots(figsize=(20,20)) 
forviolins={}
for i in range(len(set(normalizedadata.obs[key1]))):
    if '65' in list(set(normalizedadata.obs[key1]))[i]:
        cellinds=['Astrocyte' in c for c in set(normalizedadata.obs[key2])]
        forviolins[list(set(normalizedadata.obs[key1]))[i]]=varmat[i,cellinds,:].flatten()
sns.violinplot(data=pd.DataFrame(forviolins),inner='box')
plt.xticks(rotation=45)
plt.savefig(os.path.join(sc.settings.figdir,key1+key2+'AstrocyteViolins.pdf'))
plt.close()

fig, ax = plt.subplots(figsize=(20,20)) 
forviolins={}
for i in range(len(set(normalizedadata.obs[key1]))):
    if '65' in list(set(normalizedadata.obs[key1]))[i]:
        cellinds=['Astrocyte' in c for c in set(normalizedadata.obs[key2])]
        forviolins[list(set(normalizedadata.obs[key1]))[i]]=meanmat[i,cellinds,:].flatten()
sns.violinplot(data=pd.DataFrame(forviolins),inner='box')
plt.xticks(rotation=45)
plt.savefig(os.path.join(sc.settings.figdir,key1+key2+'AstrocyteVarindexViolins.pdf'))
plt.close()


fig, ax = plt.subplots(figsize=(20,20)) 
forviolins={}
for i in range(len(set(normalizedadata.obs[key1]))):
    if '65' in list(set(normalizedadata.obs[key1]))[i]:
        cellinds=['Neuron' in c for c in set(normalizedadata.obs[key2])]
        forviolins[list(set(normalizedadata.obs[key1]))[i]]=varmat[i,cellinds,:].flatten()
sns.violinplot(data=pd.DataFrame(forviolins),inner='box')
plt.xticks(rotation=45)
plt.savefig(os.path.join(sc.settings.figdir,key1+key2+'NeuronViolins.pdf'))
plt.close()

fig, ax = plt.subplots(figsize=(20,20)) 
forviolins={}
for i in range(len(set(normalizedadata.obs[key1]))):
    if '65' in list(set(normalizedadata.obs[key1]))[i]:
        cellinds=['Neuron' in c for c in set(normalizedadata.obs[key2])]
        forviolins[list(set(normalizedadata.obs[key1]))[i]]=meanmat[i,cellinds,:].flatten()
sns.violinplot(data=pd.DataFrame(forviolins),inner='box')
plt.xticks(rotation=45)
plt.savefig(os.path.join(sc.settings.figdir,key1+key2+'NeuronVarindexViolins.pdf'))
plt.close()

fig, ax = plt.subplots(figsize=(20,20)) 
forviolins={}
for i in range(len(set(normalizedadata.obs[key1]))):
    if '65' in list(set(normalizedadata.obs[key1]))[i]:
        cellinds=['Microglia' in c for c in set(normalizedadata.obs[key2])]
        forviolins[list(set(normalizedadata.obs[key1]))[i]]=varmat[i,cellinds,:].flatten()
sns.violinplot(data=pd.DataFrame(forviolins),inner='box')
plt.xticks(rotation=45)
plt.savefig(os.path.join(sc.settings.figdir,key1+key2+'MicrogliaViolins.pdf'))
plt.close()


for r in set(normalizedadata.obs['classname']):
    a=normalizedadata[normalizedadata.obs['classname']==r,:].copy()
    sc.pp.highly_variable_genes(a, n_top_genes=4000)
    sc.pp.scale(a, max_value=10)
    sc.pp.regress_out(adata,'batch')
    sc.pp.pca(a)
    sc.pp.neighbors(a)
    sc.tl.umap(a)
    sc.tl.leiden(a)
    sc_utils.cell_cycle_score(a)
    sc.pl.umap(a, color=['leiden','phase'],save='_'+r+'_summary')
    sc.pl.umap(a, color=['region'],save='_'+r+'_region')
    sc.pl.umap(a, color=['tp'],save='_'+r+'_tp')
    sc.pl.umap(a, color=['ActualRegion'],save='_'+r+'_ActualRegion')
    a.obs=adata.obs.loc[a.obs.index,:]
    sc.pl.umap(a, color=['subclassname'],save='_'+r+'_Subclass')

###Now the velocity, using the already projected data

velodirs=["E65-2019A_AND_E65-2019B_MULTI-SEQ_1_Out_velocyto","E65-2019A_AND_E65-2019B_MULTI-SEQ_2_Out_velocyto","E65-2019A_AND_E65-2019B_MULTI-SEQ_3_Out_velocyto","E80-2019_MULTI-SEQ_Out_velocyto","E90-2019_MULTI-SEQ_Out_velocyto"]
velofiles=[os.listdir(os.path.join(headpath,x))[0] for x in velodirs]
velopaths=[os.path.join(headpath,x,y) for x,y in zip(velodirs,velofiles)]
velolooms=[scv.read(x,cache=True) for x in velopaths]
for i in range(len(velolooms)):
    velolooms[i].obs.index=[re.sub("-1","",x) for x in velolooms[i].obs.index]
    velolooms[i].obs.insert(0,'batch',velodirs[i])
vdata=sc.AnnData.concatenate(*velolooms)
print(vdata.obs)
vdata.var_names_make_unique()
avdata=scv.utils.merge(adata,vdata)
print(avdata.obs)
print(adata.obs)

avdata.var_names_make_unique()
print('norm')
scv.pp.filter_genes(avdata)
scv.pp.normalize_per_cell(avdata)
scv.pp.filter_genes_dispersion(avdata)
scv.pp.log1p(avdata)
print(avdata)
print('moment')
scv.pp.moments(avdata, n_pcs=30, n_neighbors=30)
print('velo')
scv.tl.umap(avdata)
scv.tl.velocity(avdata)
print('graph')
scv.tl.velocity_graph(avdata)
scv.tl.velocity_embedding(avdata, basis='umap')
scv.pl.velocity_embedding(avdata, basis='umap',save='_embed')
scv.pl.velocity_embedding_grid(avdata, basis='umap',save='_grid')
scv.pl.velocity_embedding_stream(avdata, basis='umap',save='stream')
sc.tl.leiden(avdata)
#sc.pl.umap(avdata, color=['leiden','tp','batch'],save='_tp')
#sc.pl.umap(avdata, color=['ActualRegion'],save='_ActualRegion')
#sc.pl.umap(avdata, color=['region'],save='_region')



sc.tl.rank_genes_groups(adata, 'leiden', method='logreg')
result = adata.uns['rank_genes_groups']
groups = result['names'].dtype.names
df=pd.DataFrame({group + '_' + key[:1]: result[key][group] for group in groups for key in ['names', 'scores']})
df.to_csv(os.path.join(sc.settings.figdir,"leidenLogRegMarkers.csv"))

sc.tl.rank_genes_groups(adata, 'ActualRegion', method='logreg')
result = adata.uns['rank_genes_groups']
groups = result['names'].dtype.names
df=pd.DataFrame({group + '_' + key[:1]: result[key][group] for group in groups for key in ['names', 'scores']})
df.to_csv(os.path.join(sc.settings.figdir,"RegionLogRegMarkers.csv"))

sc.tl.rank_genes_groups(adata, 'tp', method='logreg')
result = adata.uns['rank_genes_groups']
groups = result['names'].dtype.names
df=pd.DataFrame({group + '_' + key[:1]: result[key][group] for group in groups for key in ['names', 'scores']})
df.to_csv(os.path.join(sc.settings.figdir,"tpLogRegMarkers.csv"))

