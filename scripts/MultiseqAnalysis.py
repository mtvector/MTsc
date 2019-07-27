#!/usr/bin/env python
# coding: utf-8

# In[1]:

import os
import scanpy
import scanpy as sc
sc.settings.verbosity = 3  # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.settings.set_figure_params(dpi=300)  # dots (pixels) per inch determine size of inline figures
sc.settings.cachedir=os.path.expanduser('~/macaquecache/')
sc.logging.print_versions()
import pandas as pd
import numpy as np
import re
from collections import Counter
import seaborn
import matplotlib.pyplot as plt
import anndata


# In[2]:
mspath=os.path.expanduser('~/MULTISEQmacaque')
sc.settings.figdir=os.path.expanduser('~/figs/')
sc.settings.autosave=True
sc.settings.autoshow=False

mskey=pd.read_csv(mspath+'/MultiseqKey.txt',sep='\t')
fullname=[]
mskey=mskey.astype(str,copy=True)
for i in range(mskey.shape[0]):
    if "nan" not in mskey.iloc[i,2] :
        fullname.append(mskey.iloc[i,1]+"_"+mskey.iloc[i,2]+"_"+mskey.iloc[i,3])
    else:
        fullname.append(mskey.iloc[i,1]+"_"+str(mskey.iloc[i,3]))
mskey['fullname']=fullname


# In[3]:


mskey


# In[4]:


samples=["E65-1","E65-2","E65-3","E80","E90"]
sampletable={}
for s in samples:
    sampletable[s]=pd.read_csv(mspath+'/MULTIseqOutputs/pollena-MULTIseq_'+s+'_finalcalls.txt',sep='\t').T

for s in sampletable.keys():
    sampletable[s].iloc[:,0]=[re.sub("Bar","",x) for x in sampletable[s].iloc[:,0]]

IndConversion=pd.read_csv(os.path.expanduser(mspath+'/8-12IndexConversion.txt'),sep='\t')
IndConversion=IndConversion.astype(str,copy=True)
ICdict=dict(zip(IndConversion["8-Index"],IndConversion["12-Index"]))
ICdict['Doublet']='Doublet'
ICdict['Negative']='Negative'

mskey.iloc[:,4]=[ICdict[x] for x in mskey.iloc[:,4]]

#for s in sampletable.keys():
#    sampletable[s].iloc[:,0]=[ICdict[x] for x in sampletable[s].iloc[:,0]]


# In[5]:


mskey


# In[6]:


for s in sampletable.keys():
    seaborn.countplot(sampletable[s][1]).set_title('s')
    print(sorted(Counter(sampletable[s][1]).keys()))


# In[7]:


sampleregions={}
for s in sampletable.keys():
    curkey=mskey.iloc[[x in s for x in mskey.iloc[:,1]],:]
    curkeydict=dict(zip(curkey.iloc[:,4],curkey.iloc[:,6]))
    print(curkey)
    samplematched=[]
    for i in sampletable[s].iloc[:,0]:
        inds=np.where([x == i for x in curkey.iloc[:,4]])[0]
        if len(inds)<1:
            samplematched.append(i)
        else:
            samplematched.append(curkey.iat[ inds[0] ,6])
    sampletable[s].iloc[:,0]=samplematched


# In[8]:


for s in sampletable.keys():
    seaborn.countplot(sampletable[s][1]).set_title('s')
    print(sorted(Counter(sampletable[s][1]).keys()))


# In[10]:


def cell_cycle_score(adata,save=False):
    s_genes=['MCM5', 'PCNA', 'TYMS', 'FEN1', 'MCM2', 'MCM4', 'RRM1', 'UNG', 'GINS2', 'MCM6', 'CDCA7', 'DTL', 'PRIM1', 'UHRF1', 'MLF1IP', 'HELLS', 'RFC2', 'RPA2', 'NASP', 'RAD51AP1', 'GMNN', 'WDR76', 'SLBP', 'CCNE2', 'UBR7', 'POLD3', 'MSH2', 'ATAD2', 'RAD51', 'RRM2', 'CDC45', 'CDC6', 'EXO1', 'TIPIN', 'DSCC1', 'BLM', 'CASP8AP2', 'USP1', 'CLSPN', 'POLA1', 'CHAF1B', 'BRIP1', 'E2F8']
    g2m_genes=['HMGB2', 'CDK1', 'NUSAP1', 'UBE2C', 'BIRC5', 'TPX2', 'TOP2A', 'NDC80', 'CKS2', 'NUF2', 'CKS1B', 'MKI67', 'TMPO', 'CENPF', 'TACC3', 'FAM64A', 'SMC4', 'CCNB2', 'CKAP2L', 'CKAP2', 'AURKB', 'BUB1', 'KIF11', 'ANP32E', 'TUBB4B', 'GTSE1', 'KIF20B', 'HJURP', 'CDCA3', 'HN1', 'CDC20', 'TTK', 'CDC25C', 'KIF2C', 'RANGAP1', 'NCAPD2', 'DLGAP5', 'CDCA2', 'CDCA8', 'ECT2', 'KIF23', 'HMMR', 'AURKA', 'PSRC1', 'ANLN', 'LBR', 'CKAP5', 'CENPE', 'CTCF', 'NEK2', 'G2E3', 'GAS2L3', 'CBX5', 'CENPA']
    sc.tl.score_genes_cell_cycle(adata,s_genes=s_genes,g2m_genes=g2m_genes)
    sc.pl.violin(adata, ['G2M_score','S_score'], groupby='leiden')
    if 'X_tsne' in adata.obsm.keys(): 
        sc.pl.tsne(adata, color=['G2M_score','S_score','phase','leiden'],save='_cc')
    sc.pl.umap(adata, color=['G2M_score','S_score','phase','leiden'],save='_cc')
    return(adata)

# In[18]:


sampleDirs=dict(zip(sampletable.keys(),['E65-2019A_AND_E65-2019B_MULTI-SEQ_1_Out',
                                        'E65-2019A_AND_E65-2019B_MULTI-SEQ_2_Out',
                                        'E65-2019A_AND_E65-2019B_MULTI-SEQ_3_Out',
                                        'E80-2019_MULTI-SEQ_Out',
                                       'E90-2019_MULTI-SEQ_Out']))


# In[19]:


adatas=[]
for s in sampletable.keys():
    adata = sc.read_10x_mtx('/scrapp2/mtschmitz/'+sampleDirs[s]+'/outs/filtered_feature_bc_matrix',cache=True)
    adata.obs.index=[re.sub("-1","",x) for x in adata.obs.index]
    sampletable[s].columns=['region']
    adata.obs.insert(0,'region',sampletable[s].loc[adata.obs.index,'region'])
    adatas.append(adata)
    


# In[20]:


adata=anndata.AnnData.concatenate(*adatas)

def tp_format_macaque(name):
    if isinstance(name, str):
        searched=re.search("E[0-9]+",name)
        if searched is None:
            return('nan')
        tp=re.sub("^E","",searched.group(0))
        tp=float(tp)
        return(tp)
    else:
        return('nan')

def region_format_macaque(name):
    if isinstance(name, str):
        region=re.sub('E[0-9]+',"",name)
        region=re.sub("^_","",region)
        if '_' not in name:
            return('nan')
        spl=region.split('_')
        if len(spl) < 1:
            return('nan')
        region=spl[len(spl)-1]
    else:
        region='nan'
    return(region)


# In[25]:

adata.var_names_make_unique()

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
sc.pl.violin(adata, ['n_genes', 'n_counts', 'percent_mito'],
             jitter=0.4, multi_panel=True)

adata.obs.insert(0, 'simpleregion',[re.sub('_','',x[x.rfind('_'):]) if '_' in x else x for x in adata.obs['region']])
adata.obs.insert(0, 'tp',[tp_format_macaque(x) for x in adata.obs['region']])
adata.obs.insert(0, 'ActualRegion',[region_format_macaque(x) for x in adata.obs['region']])

sc.pp.normalize_per_cell(adata, counts_per_cell_after=1e4)

sc.pp.log1p(adata)

sc.pp.highly_variable_genes(adata, n_top_genes=10000)
sc.pp.scale(adata, max_value=10)
#sc.pp.regress_out(adata,'n_counts')

sc.pp.pca(adata)
sc.pp.neighbors(adata)
sc.tl.umap(adata)
sc.tl.leiden(adata)


#cell_cycle_score(adata)
sc.pl.umap(adata, color=['leiden','n_counts','region'])
sc.pl.umap(adata, color=['region','simpleregion','batch','leiden'],save='_MultiseqSummary')
sc.pl.umap(adata, color=['region'],save='_region')
sc.pl.umap(adata, color=['leiden'],save='_leiden')
sc.pl.umap(adata, color=['tp'],save='_tp')
sc.pl.umap(adata, color=['ActualRegion'],save='_ActualRegion')


sc.pl.umap(adata[['Negative' not in x and 'Doublet' not in x for x in adata.obs['region']],:], color=['leiden'],save='_screenedleiden')
sc.pl.umap(adata[['Negative' not in x and 'Doublet' not in x for x in adata.obs['region']],:], color=['region'],save='_screenedsample')
sc.pl.umap(adata[['Negative' not in x and 'Doublet' not in x for x in adata.obs['region']],:], color=['simpleregion'],save='_screenedsimpleregion')


# In[26]:


#sc.tl.rank_genes_groups(adata, 'leiden', method='logreg')
#result = adata.uns['rank_genes_groups']
#groups = result['names'].dtype.names
#df=pd.DataFrame({group + '_' + key[:1]: result[key][group] for group in groups for key in ['names', 'scores']})
#df.to_csv(os.path.join(sc.settings.figdir,"leidenLogRegMarkers.csv"))

def marker_analysis(adata,variables=['leiden','region'],markerpath='https://docs.google.com/spreadsheets/d/e/2PACX-1vTz5a6QncpOOO-f3FHW2Edomn7YM5mOJu4z_y07OE3Q4TzcRr14iZuVyXWHv8rQuejzhhPlEBBH1y0V/pub?gid=1154528422&single=true&output=tsv'):
    sc.set_figure_params(color_map="Purples")
    import random
    markerpath=os.path.expanduser(markerpath)
    markers=pd.read_csv(markerpath,sep="\t")
    print(markers)
    markers[markers.keys()[0]]=[str(x) for x in markers[markers.keys()[0]]]
    markers[markers.keys()[2]]=[str(x).split(',') for x in markers[markers.keys()[2]]]
    markers[markers.keys()[3]]=[str(x).split(';') for x in markers[markers.keys()[3]]]
    markers[markers.keys()[3]]=[[str(x).split(',') for x in y] for y in markers[markers.keys()[3]]]
    uniqueClasses=set([y for x in markers[markers.keys()[2]] for y in x if y!='nan'])
    uniqueSubClasses=set([z for x in markers[markers.keys()[3]] for y in x for z in y if z!='nan'])
    comboClasses=[]
    print(markers)
    for i in range(markers.shape[0]):
        rowlist=[]
        for j in range(len(markers[markers.keys()[2]][i])):
            for k in markers[markers.keys()[3]][i][j]:
                rowlist.append(' '.join(filter(lambda x: x != 'nan',[k,markers[markers.keys()[2]][i][j]])))
        comboClasses.append(rowlist)
    markers['fullclass']=comboClasses
    markers.set_index(markers.keys()[0],inplace=True,drop=False)
    markers=markers.loc[ [x for x in markers[markers.keys()[0]] if x in adata.var_names],:]
    uniqueFullClasses=set([y for x in markers['fullclass'] for y in x if y!='nan'])
    from collections import defaultdict
    markerDict = defaultdict(list)

    for x in uniqueFullClasses:
        for y in markers[markers.keys()[0]]:
            if x in markers.loc[y,'fullclass']:
                markerDict[x].append(y)
    markerDictClass = defaultdict(list)
    for x in uniqueClasses:
        for y in markers[markers.keys()[0]]:
            if x in markers.loc[y,'fullclass']:
                markerDictClass[x].append(y)

    markerPlotGroups=[]
    for k in markerDict.keys():
        if len(markerDict[k])>1:
            print(k)
            print(len(markerDict[k]))
            sc.tl.score_genes(adata,gene_list=markerDict[k],score_name=k,gene_pool= markerDict[k]+random.sample(adata.var.index.tolist(),min(4000,adata.var.index.shape[0])))
            markerPlotGroups.append(k)
    adata.uns['marker_groups']=list(markerDict.keys())
    for tag in variables:
        pd.DataFrame(adata.obs.groupby(tag).describe()).to_csv(os.path.join(sc.settings.figdir, tag+"MarkerSumStats.csv"))

    if 'X_tsne' in adata.obsm.keys():
        sc.pl.tsne(adata, color=markerPlotGroups,save="_Marker_Group")
    sc.pl.umap(adata, color=markerPlotGroups,save="_Marker_Group")
    print(markerDict)
    #sc.pl.violin(adata, markerPlotGroups, groupby='leiden',save="_Marker_Group_violins")
    for i in markerDictClass:
        print(i)
        if 'X_tsne' in adata.obsm.keys():
            sc.pl.tsne(adata, color=sorted(markerDictClass[i]),save="_"+str(i)+"_Marker")
        sc.pl.umap(adata, color=sorted(markerDictClass[i]),save="_"+str(i)+"_Marker")
    return(adata)

marker_analysis(adata,variables=['leiden','region'],markerpath=os.path.expanduser('~/markers.txt'))

import scvelo as scv
scv.settings.figdir=os.path.expanduser('~/figs/')
scv.settings.autosave=True
scv.settings.autoshow=False

headpath=os.path.expanduser("/scrapp2/mtschmitz")
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
scv.pp.filter_and_normalize(avdata)
scv.tl.velocity_graph(avdata)
scv.tl.velocity_embedding(avdata, basis='umap')
scv.pl.velocity_embedding(avdata, basis='umap',save='Embed')
scv.pl.velocity_embedding_grid(avdata, basis='umap',save='Grid')
