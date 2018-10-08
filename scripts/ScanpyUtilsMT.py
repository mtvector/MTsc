import os
import scanpy
import scanpy.api as sc
import pandas as pd
import numpy as np
import sklearn
from sklearn import decomposition
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
matplotlib.use('Agg')
plt.switch_backend('agg')
import logging as logg
import seaborn
import pickle
import inspect
import re

sc.set_figure_params(color_map="copper")
outpath=sc.settings.figdir
#Init scanpyanalysis object by giving an outpath, then call functions statically
outpath = os.path.expanduser(outpath)
sc.settings.autosave=True
sc.settings.autoshow=False
if not os.path.exists(outpath):
    os.makedirs(outpath)

#Saves a scanpy object from this pipeline (includes name and list of operations run in uns) as a h5ad
def save_adata(adata,outpath=None):
    if outpath is None:
        outpath=sc.settings.figdir
    adata.write(os.path.join(outpath,"-".join([adata.uns['name']]+list(adata.uns['operations']))+".h5ad"))

#Template function for pipeline
def template(adata,save=False):
    #anndata.read_h5ad(file)
    adata.uns['operations']=np.append(adata.uns['operations'],inspect.stack()[0][3])
    if save:
        save_adata(adata,outpath)
    return(None)

#import a single 10x file, do basic
def import_file(filename,refname,outpath=None,filter_genes=30,save=False):
    if outpath is None:
        outpath=sc.settings.figdir
    adata=sc.read_10x_h5(filename,refname)
    adata.var_names_make_unique()
    #sc.pp.filter_cells(adata, min_genes=filter_genes)
    adata.obs['n_counts'] = adata.X.sum(axis=1).A1
    print(adata)
    adata.uns['name']=filename.split(os.sep)[-3]
    adata.uns['operations']=['load']
    if save:
        save_adata(adata,outpath)
    return(adata)

#Use these functions to modularize munging of tp, region names from name string
#These can be remade for other applications, or generalized so you can pass a list of functions
def tp_format_macaque(name):
    return(float(re.sub("^E","",re.search("E[0-9]+",name).group(0))))

def region_format_macaque(name):
    region=re.sub('E[0-9]+',"",name)
    region=re.sub("^_","",region)
    region=region.split('_')[0]
    return(region)

def macaque_process_irregular_names(adata):
    mixSamples={'Mix1_Out':"E100caudateANDextBGANDintBGANDputanum",
    'Mix2_Out':"E100pulvinarANDdorsalThalamusANDventralThalamus",
    'Mix3_Out':"E100lateralANDanteriorANDvermis",
    'Mix4_Out':"E80lateralANDanteriorANDvermis",
    'Mix5_Out':"E80putanumANDclaustrumANDbasalGanglia",
    'Mix6_Out':"E80dorsalThalamusANDventralThalamus",
    'Mix7_Out':"E100hypothalamusANDPoA",
    'Mix8_Out':"E100septumANDnuclearAccumbens",
    'Mix9_Out':"E100CGE-LGE",
    'Mix10_Out':"E80CGE-LGE",
    'Mix11_Out':"E90choroid",
    'E10somato_Out':'E100somato',
    'E100insulaCP_Out':'E100insula',
    'E40_hippocampus_Out':'E40hippo',
    'E50_hippo-temp_Out':'E50hippo',
    'E65_hippocampus_Out':'E65hippo',
    'E65_somatosensory_Out':'E65somato',
    'E65_hypothalamus_Out':'E65hypo'}
    if adata.uns['name'] in mixSamples.keys():
        adata.uns['name']=mixSamples[adata.uns['name']]
    return(adata)

#Imports a list of files, filters them with rough filtering, returns list of scanpy objects.
#For use with concat_files function
#Set filters functions to None if not processing developing macaque brain
def import_files(fileList,refname,groupname="group",n_counts=400,percent_mito=.6,percent_ribo=.6,filter_genes=10,tp_func=tp_format_macaque,region_func=region_format_macaque,fix_irreg=macaque_process_irregular_names,log=False,save=True):
    adatas=[]
    batch_categories=[]
    for filename in fileList:
        a = sc.read_10x_h5(filename,refname)
        a.uns['operations']=['load_file']
        sc.pp.filter_cells(a, min_genes=filter_genes)
        a.obs['n_counts'] = a.X.sum(axis=1).A1
        a=quantify_ribo(a)
        a=quantify_mito(a)
        a=filter_custom(a,n_counts=n_counts, percent_mito=percent_mito,percent_ribo=percent_ribo)
        #a=knee_filter_cells(a)
        if log:
            sc.pp.log1p(a)
        a.uns['name']=filename.split(os.sep)[-3]
        a.obs['batch']=a.uns['name']
        if fix_irreg is not None:
            a=fix_irreg(a)
        if tp_func is not None:
            print(a)
            print(a.uns['name'])
            a.obs['timepoint']=tp_func(a.uns['name'])
        if region_func is not None:
            a.obs['region']=region_func(a.uns['name'])
        batch_categories.append(a.uns['name'])
        adatas.append(a)
    return(adatas)

#merges a list of scanpy objects into a single object
#This is a starting point for the pipeline
def concat_files(fileList,refname,outpath=None,groupname="group",n_counts=400,percent_mito=.6,percent_ribo=.6,filter_genes=10,save=True):
    if outpath is None:
        outpath=sc.settings.figdir
    adatas=import_files(fileList=fileList,refname=refname,groupname=groupname,n_counts=n_counts,percent_mito=percent_mito,percent_ribo=percent_ribo,filter_genes=filter_genes)
    adata = sc.AnnData.concatenate(*adatas)
    adata.uns['name']=groupname
    adata.uns['operations']=['concat_files']
    sc.settings.autosave=True
    sc.settings.autoshow=False
    print(adata)
    if save:
        save_adata(adata,outpath)
    return(adata)

#Not working yet
def concat_files_mnn(fileList,refname,outpath=None,groupname="group",n_top_genes=None,n_counts=400,percent_mito=.6,percent_ribo=.6,filter_genes=10,save=True):
    if outpath is None:
        outpath=sc.settings.figdir
    if n_top_genes is None:
        n_top_genes=np.sum(np.sum(adata.X,axis=1)>0)
    adata=concat_files(fileList=fileList,refname=refname,groupname="macaqueDevBrain",n_counts=500,percent_mito=.5,percent_ribo=.5,filter_genes=10,save=True)
    var_subset=sc.pp.filter_genes_dispersion(  # select highly-variable genes
        adata.X, flavor='seurat', n_top_genes=n_top_genes, log=True)
    var_subset=adata.var.index[var_subset.gene_subset]
    print(adata)
    print(adata.X)
    print(var_subset)
    adatas=import_files(fileList=fileList,refname=refname,groupname=groupname,n_counts=n_counts,percent_mito=percent_mito,percent_ribo=percent_ribo,filter_genes=filter_genes,log=True)
    adata = sc.pp.mnn_correct(*adatas,var_index=var_subset)
    adata.uns['name']=groupname+"_MNN"
    adata.uns['operations']=np.append(adata.uns['operations'],inspect.stack()[0][3])
    if save:
        save_adata(adata,outpath)
    return(adata)

#Basic filtering functions
def filter_custom(adata,percent_mito=.5,percent_ribo=.5,n_counts=400):
    #adata.write(os.path.expanduser("~/Desktop/tmp.h5ad"))
    #adata=adata[adata.obs.n_counts>n_counts,:]
    adata._inplace_subset_obs(adata.obs.n_counts>n_counts)
    #adata=adata[adata.obs.percent_mito<percent_mito,:]
    adata._inplace_subset_obs(adata.obs.percent_mito<percent_mito)
    #adata=adata[adata.obs.percent_ribo<percent_ribo,:]
    adata._inplace_subset_obs(adata.obs.percent_ribo<percent_ribo)
    adata.raw=adata
    print(adata)
    return(adata)

#Count the percent of ribosomal genes per cell, add to obs
def quantify_ribo(adata,ribo_genes=None,save=False):
    if ribo_genes is None:
        ribo_genes=[name for name in adata.var_names if name.startswith('RPS') or name.startswith('RPL') ]
    adata.obs['percent_ribo'] = np.sum(
        adata[:, ribo_genes].X, axis=1).A1 / np.sum(adata.X, axis=1).A1
    adata.uns['operations']=np.append(adata.uns['operations'],inspect.stack()[0][3])
    if save:
        save_adata(adata,sc.settings.figdir)
    return(adata)

#Count the percent of mitochondrial genes per cell, add to obs
def quantify_mito(adata,mito_genes=None,save=False):
    if mito_genes is None:
        mito_genes = [name for name in adata.var_names if name in ['ND1','ND2','ND4L','ND4','ND5','ND6','ATP6','ATP8','CYTB','COX1','COX2','COX3']]
    adata.obs['percent_mito'] = np.sum(
        adata[:, mito_genes].X, axis=1).A1 / np.sum(adata.X, axis=1).A1
    adata.uns['operations']=np.append(adata.uns['operations'],inspect.stack()[0][3])
    if save:
        save_adata(adata,sc.settings.figdir)
    return(adata)

#Doesn't work properly on macaque data
def knee_filter_cells(adata,expect_cells=30000,save=False):
    knee_cells=list(getKneeEstimate(Counter(adata.obs.n_counts.to_dict()),expect_cells=expect_cells,plotfile_prefix=sc.settings.figdir))
    inthresh=np.array([x in knee_cells for x in adata.obs.index])
    if save:
        save_adata(adata,sc.settings.figdir)
    return(adata._inplace_subset_obs(inthresh))

#Standard normalizing and scaling, used before score_genes and PCA
def norm_and_scale(adata,n_top_genes=None,min_cells=5,save=False):
    sc.pp.filter_genes(adata, min_counts=1)
    sc.pp.filter_genes(adata, min_cells=min_cells)
    if n_top_genes is None:
        n_top_genes=adata.shape[1]
    filter_result = sc.pp.filter_genes_dispersion( adata.X, flavor='cell_ranger', n_top_genes=n_top_genes, log=False)
    adata = adata[:, filter_result.gene_subset]     # subset the genes
    sc.pp.normalize_per_cell(adata)
    sc.pp.log1p(adata)
    sc.pp.scale(adata)
    adata.uns['operations']=np.append(adata.uns['operations'],inspect.stack()[0][3])
    if save:
        save_adata(adata,sc.settings.figdir)
    return(adata)

#Runs scanpy pca
def do_pca(adata,save=False):
    sc.tl.pca(adata)
    adata.obsm['X_pca'] *= -1  # multiply by -1 to match Seurat
    sc.pl.pca_scatter(adata)
    sc.pl.pca_variance_ratio(adata, log=True)
    sc.pl.pca_loadings(adata)
    adata.uns['operations']=np.append(adata.uns['operations'],inspect.stack()[0][3])
    if save:
        save_adata(adata,sc.settings.figdir)
    return(adata)

#Runs tsne, use after quantifying ribo and mito
def do_tsne(adata,save=False):
    sc.tl.tsne(adata, random_state=2, n_pcs=50)
    sc.pl.tsne(adata, color=['n_counts','percent_ribo','percent_mito'],  save="_percentRiboMito")
    adata.uns['operations']=np.append(adata.uns['operations'],inspect.stack()[0][3])
    if save:
        save_adata(adata,sc.settings.figdir)
    return(adata)

#Runs neighbors to build KNN graph from PCA, also runs UMAP
def neighbors_and_umap(adata,n_neighbors=100,save=False):
    sc.pp.neighbors(adata, n_neighbors=n_neighbors)
    sc.tl.umap(adata)
    adata.uns['operations']=np.append(adata.uns['operations'],inspect.stack()[0][3])
    if save:
        save_adata(adata,sc.settings.figdir)
    return(adata)

#Runs louvain clustering
def louvain(adata,save=False):
    sc.tl.louvain(adata)
    adata.uns['operations']=np.append(adata.uns['operations'],inspect.stack()[0][3])
    sc.pl.umap(adata,color="louvain",save="_louvain")
    if 'X_tsne' in adata.obsm.keys():
        sc.pl.tsne(adata,color="louvain",  save="_louvain")
    else:
        sc.pl.umap(adata,color="louvain",  save="_louvain")
    if save:
        save_adata(adata,sc.settings.figdir)
    return(adata)

#Function wraps together standard basic processing
def std_norm_transform(adata,n_top_genes=6000,log=True,tsne=True,save=False):
    if log:
        sc.pp.recipe_zheng17(adata,n_top_genes=n_top_genes)
    else:
         sc.pp.recipe_zheng17(adata,n_top_genes=n_top_genes,log=False)
    print(adata)
    #adata=norm_and_scale(adata,n_top_genes)
    adata=do_pca(adata)
    if tsne:
        adata=do_tsne(adata)
    adata=neighbors_and_umap(adata)
    adata=louvain(adata)
    if save:
        save_adata(adata,sc.settings.figdir)
    return(adata)

#To Do, incomplete
def woublet(adata,save=False):
    print(adata.shape[0]/100000,"doublets")
    sc.tl.woublet(adata,expected_doublet_rate=adata.shape[0]/(100000))
    adata.uns['operations']=np.append(adata.uns['operations'],inspect.stack()[0][3])
    if save:
        save_adata(adata,sc.settings.figdir)
    return(adata)

#Strips genes from the expression matrix of a scanpy objects
def strip_genes(adata,genes,groupname='genegroup',save=False):
    genes=np.array([x for x in genes if x in adata.var.index])
    print(genes)
    print(np.array(adata[:,genes].X))
    adata.obsm[str(groupname)]=np.array(adata[:,genes].X.todense())
    adata.uns[groupname+'_genenames']=genes
    adata=adata[:,[x for x in adata.var.index if x not in genes]]
    adata.uns['operations']=np.append(adata.uns['operations'],inspect.stack()[0][3])
    if save:
        save_adata(adata,sc.settings.figdir)
    return(adata)

#Calculates the cell cycle score using Aviv Regevs cell cycle stage markers
def cell_cycle_score(adata,save=False):
    s_genes=['MCM5', 'PCNA', 'TYMS', 'FEN1', 'MCM2', 'MCM4', 'RRM1', 'UNG', 'GINS2', 'MCM6', 'CDCA7', 'DTL', 'PRIM1', 'UHRF1', 'MLF1IP', 'HELLS', 'RFC2', 'RPA2', 'NASP', 'RAD51AP1', 'GMNN', 'WDR76', 'SLBP', 'CCNE2', 'UBR7', 'POLD3', 'MSH2', 'ATAD2', 'RAD51', 'RRM2', 'CDC45', 'CDC6', 'EXO1', 'TIPIN', 'DSCC1', 'BLM', 'CASP8AP2', 'USP1', 'CLSPN', 'POLA1', 'CHAF1B', 'BRIP1', 'E2F8']
    g2m_genes=['HMGB2', 'CDK1', 'NUSAP1', 'UBE2C', 'BIRC5', 'TPX2', 'TOP2A', 'NDC80', 'CKS2', 'NUF2', 'CKS1B', 'MKI67', 'TMPO', 'CENPF', 'TACC3', 'FAM64A', 'SMC4', 'CCNB2', 'CKAP2L', 'CKAP2', 'AURKB', 'BUB1', 'KIF11', 'ANP32E', 'TUBB4B', 'GTSE1', 'KIF20B', 'HJURP', 'CDCA3', 'HN1', 'CDC20', 'TTK', 'CDC25C', 'KIF2C', 'RANGAP1', 'NCAPD2', 'DLGAP5', 'CDCA2', 'CDCA8', 'ECT2', 'KIF23', 'HMMR', 'AURKA', 'PSRC1', 'ANLN', 'LBR', 'CKAP5', 'CENPE', 'CTCF', 'NEK2', 'G2E3', 'GAS2L3', 'CBX5', 'CENPA']
    sc.tl.score_genes_cell_cycle(adata,s_genes=s_genes,g2m_genes=g2m_genes)
    sc.pl.violin(adata, ['G2M_score','S_score'], groupby='louvain')
    if 'X_tsne' in adata.obsm.keys():
        sc.pl.tsne(adata, color=['G2M_score','S_score','phase','louvain'])
    else:
        sc.pl.umap(adata, color=['G2M_score','S_score','phase','louvain'])
    adata.uns['operations']=np.append(adata.uns['operations'],inspect.stack()[0][3])
    if save:
        save_adata(adata,sc.settings.figdir)
    return(adata)

#Calculates logistic regression enrichment of louvain clusters
def log_reg_diff_exp(adata,save=False):
    sc.tl.rank_genes_groups(adata, 'louvain', method='logreg')
    result = adata.uns['rank_genes_groups']
    groups = result['names'].dtype.names
    df=pd.DataFrame({group + '_' + key[:1]: result[key][group] for group in groups for key in ['names', 'scores']})
    df.to_csv(os.path.join(sc.settings.figdir,"LouvainLogRegMarkers.csv"))
    adata.uns['operations']=np.append(adata.uns['operations'],inspect.stack()[0][3])
    if save:
        save_adata(adata,sc.settings.figdir)
    return(adata)

#Carry out brain atlas marker expresion, preprocessing if very complicated to process class and subclasses of markers
#Uses score_genes
def marker_analysis(adata,markerpath='https://docs.google.com/spreadsheets/d/e/2PACX-1vTz5a6QncpOOO-f3FHW2Edomn7YM5mOJu4z_y07OE3Q4TzcRr14iZuVyXWHv8rQuejzhhPlEBBH1y0V/pub?gid=1154528422&single=true&output=tsv',save=False):
    sc.set_figure_params(color_map="Purples")
    import random
    markerpath=os.path.expanduser(markerpath)
    markers=pd.read_csv(markerpath,sep="\t")
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

    pd.DataFrame(adata.obs.groupby(['louvain']).describe()).to_csv(os.path.join(sc.settings.figdir,"ClusterMarkerSumStats.csv"))
    #This belongs outside of pipeline function
    pd.DataFrame(adata.obs.groupby(['region']).describe()).to_csv(os.path.join(sc.settings.figdir,"RegionMarkerSumStats.csv"))

    if 'X_tsne' in adata.obsm.keys():
        sc.pl.tsne(adata, color=markerPlotGroups,save="_Marker_Group")
    sc.pl.umap(adata, color=markerPlotGroups,save="_Marker_Group")

    sc.pl.violin(adata, markerPlotGroups, groupby='louvain',save="_Marker_Group_violins")
    for i in markerDictClass:
        if 'X_tsne' in adata.obsm.keys():
            sc.pl.tsne(adata, color=sorted(markerDictClass[i]),save="_"+str(i)+"_Marker")
        sc.pl.umap(adata, color=sorted(markerDictClass[i]),save="_"+str(i)+"_Marker")
    #General
    adata.uns['operations']=np.append(adata.uns['operations'],inspect.stack()[0][3])
    if save:
        save_adata(adata,sc.settings.figdir)
    return(adata)

def dirichlet_marker_analysis(adata,markerpath='https://docs.google.com/spreadsheets/d/e/2PACX-1vTz5a6QncpOOO-f3FHW2Edomn7YM5mOJu4z_y07OE3Q4TzcRr14iZuVyXWHv8rQuejzhhPlEBBH1y0V/pub?gid=1154528422&single=true&output=tsv'):
    sc.set_figure_params(color_map="Purples")
    import random
    markerpath=os.path.expanduser(markerpath)
    markers=pd.read_csv(markerpath,sep="\t")
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

    markerPosterior = pd.DataFrame()
    for k in markerDict.keys():
        if len(markerDict[k])>1:
            markerPosterior=pd.concat([markerPosterior,pd.DataFrame(np.mean(adata[:,np.array(markerDict[k])].varm['gene_topic'],axis=0),columns=[k])],axis=1)

    markerPosterior.to_csv(os.path.join(sc.settings.figdir,"MarkerPosteriorMeans.csv"))

#Calculate Hierarchical dirichlet process model for expression set
def sc_hdp(adata,alpha=1,eta=.01,gamma=1,eps=1e-5,save=False):
    import gensim
    import scipy
    import itertools
    def moving_average(data,window_width):
        cumsum_vec = np.cumsum(np.insert(data, 0, 0))
        ma_vec = (cumsum_vec[window_width:] - cumsum_vec[:-window_width]) / window_width
        return(ma_vec)

    def table_top_words(word_topic, feature_names, n_top_words):
        df=[]
        for topic_idx, topic in enumerate(word_topic):
            print(topic_idx)
            df.append([feature_names[i]
                                for i in topic.argsort()[:-n_top_words - 1:-1]])
        print(df)
        return(pd.DataFrame(df))

    #A bit of a hack because normally this is returned nice and sparse, but you don't know which topics are included
    #Works, slow, only method that works
    def get_doc_topic(corpus, model,eps=0.0):
        doc_topic = list()
        for doc in corpus:
            doc_topic.append(model.__getitem__(doc, eps=eps))
            #Convert index,value tuples to sparse matrix, then dense
        ii=[[i]*len(v) for i,v in enumerate(doc_topic)]
        ii=list(itertools.chain(*ii))
        jj=[j for j,_ in itertools.chain(*doc_topic)]
        data=[d for _,d in itertools.chain(*doc_topic)]
        return(scipy.sparse.csr_matrix((data, (ii, jj))).todense(),list(set(jj)))


    def get_topic_to_wordids(model,eps=0.0):
        p = list()
        if hasattr(model, 'm_T'):
            for topicid in range(model.m_T):
                topic = model.m_lambda[topicid]
                topic = topic / topic.sum() # normalize to probability dist
                p.append(topic)
            return(np.array(p).T)
        else:
            for topicid in range(model.num_terms):
                topic = model.get_term_topics(topicid,minimum_probability=eps)
                p.append(topic)
            ii=[[i]*len(v) for i,v in enumerate(p)]
            ii=list(itertools.chain(*ii))
            jj=[j for j,_ in itertools.chain(*p)]
            data=[d for _,d in itertools.chain(*p)]
            return(scipy.sparse.csr_matrix((data, (ii, jj))).todense(),list(set(jj)))

    #Merges 2 2D dataframes into 3d, then collapses them every other column (Name1,Val1,Name2,Val2)
    #word_topic must be np array, not np matrix
    def table_top_words(word_topic, feature_names):
        return(pd.DataFrame(np.stack([feature_names[(-word_topic).argsort(axis=0)],word_topic[(-word_topic).argsort(axis=0)][:,:,0]],2).reshape(word_topic.shape[0],-1)))

    #def table_top_words(word_topic, feature_names, n_top_words):
    #    df=[]
    #    for topic_idx, topic in enumerate(word_topic.T):
    #        print(topic_idx)
    #        df.append([feature_names[i]
    #                            for i in topic.argsort()[:-n_top_words - 1:-1]])
    #print(df)
    #return(pd.DataFrame(df))

    def intersection(lst1, lst2):
        return list(set(lst1) & set(lst2))

    adata=adata[:,adata.var.index.argsort()]
    model = gensim.models.HdpModel(corpus=gensim.matutils.Sparse2Corpus( adata.X.T),id2word=gensim.corpora.dictionary.Dictionary([adata.var.index.tolist()]),alpha=alpha,gamma=gamma,eta=eta)
    print(model)
    modelLDA=model.suggested_lda_model()
    #Will drop topics even with eps = 0.0 ... idk why
    doc_topic,included_topics=get_doc_topic(gensim.matutils.Sparse2Corpus(adata.X.T),modelLDA,eps=eps)
    print(doc_topic)
    print(np.array(doc_topic))
    print(doc_topic.shape)
    print(type(doc_topic))
    print(included_topics)
    #Included topics from words seems unnecessary empirically, but will catch in variable anyway
    word_topic,included_topics_from_words=get_topic_to_wordids(modelLDA,eps=0.0)
    doc_topic=doc_topic[:,intersection(included_topics,included_topics_from_words)]
    word_topic=word_topic[:,intersection(included_topics,included_topics_from_words)]
    #Slim down the topics to the ones that are well represented in cells (sum > eps)
    included_topics=np.array(np.sum(doc_topic,axis=0)/sum(np.sum(doc_topic,axis=0))>eps)[0]
    print(included_topics)
    print(len(included_topics))
    print(word_topic)
    print(word_topic.shape)
    doc_topic=doc_topic[:,included_topics]
    word_topic=word_topic[:,included_topics]
    adata.varm['gene_topic']=word_topic
    adata.obsm['cell_topic']=doc_topic
    adata.uns['operations']=np.append(adata.uns['operations'],inspect.stack()[0][3])
    if save:
        save_adata(adata,sc.settings.figdir)
    table_top_words(np.array(word_topic),adata.var.index).to_csv(os.path.join(sc.settings.figdir,"TopicMarkers.txt"))
    pl=seaborn.barplot(x=list(range(doc_topic.shape[1])),y=np.sum(np.array(doc_topic),axis=0)).get_figure()
    pl.savefig(os.path.join(sc.settings.figdir,"TopicWeightBarplot.png"))
    plt.clf()
    return(adata)

#Calculate non-negative matrix factorization of a dataset
def sc_nmf(adata,n_components=12,save=True):
    import sklearn
    from sklearn import decomposition
    def moving_average(data,window_width):
        cumsum_vec = np.cumsum(np.insert(data, 0, 0))
        ma_vec = (cumsum_vec[window_width:] - cumsum_vec[:-window_width]) / window_width
        return(ma_vec)
    def table_top_words(model, feature_names, n_top_words):
        df=[]
        for topic_idx, topic in enumerate(model.components_):
            print(topic_idx)
            df.append([feature_names[i]
                                for i in topic.argsort()[:-n_top_words - 1:-1]])
        print(df)
        return(pd.DataFrame(df))

    adata=adata[adata.obs.n_counts.argsort()[::-1],:]
    nmfM=sklearn.decomposition.NMF(n_components=n_components,verbose=1)
    nmfM.fit(adata.X)
    doc_topic=nmfM.transform(adata.X)
    adata.uns['nmf_model']=nmfM
    for i in range(doc_topic.shape[1]):
        adata.obs['nmf_'+str(i)]=doc_topic[:,i]

    for i in range(nmfM.components_.shape[0]):
        adata.var['nmf_'+str(i)]=np.log10(nmfM.components_[i,:])
    import seaborn
    hm=seaborn.heatmap(doc_topic).get_figure()
    hm.savefig(os.path.join(sc.settings.figdir,"TopicHeatmap.png"))
    plt.clf()
    plt.plot(range(doc_topic.shape[0]),np.std( doc_topic,axis=1))
    plt.title("SD of topic allocation in cells")
    plt.xlabel("RankedCells")
    plt.ylabel("SD")
    plt.savefig(os.path.join(sc.settings.figdir,"TopicCellSD"))
    plt.clf()
    convolvedSD=moving_average(np.std( doc_topic,axis=1),300)
    plt.plot(range(len(convolvedSD)),convolvedSD)
    plt.title("Moving Average Sparsity (SD)")
    plt.xlabel("RankedCells")
    plt.ylabel("SD")
    plt.savefig(os.path.join(sc.settings.figdir,"TopicCellMASD.png"))
    plt.clf()
    plt.plot(range(adata.obs.shape[0]),np.log10(adata.obs.n_counts))
    plt.title("Log Counts per Cell")
    plt.xlabel("RankedCells")
    plt.title("Log Reads")
    plt.savefig(os.path.join(sc.settings.figdir,"LogCounts.png"))
    plt.clf()
    for i in range(doc_topic.shape[1]):
        convolved=moving_average(doc_topic[:,i],300)
        plt.plot(range(len(convolved)),convolved)
    plt.xlabel("RankedCells")
    plt.ylabel("Moving Avg Topic Probability")
    plt.savefig(os.path.join(sc.settings.figdir,"TopicMovingAvg.png"))
    plt.clf()
    table_top_words(nmfM,adata.var.index,40).to_csv(os.path.join(sc.settings.figdir,"TopicMarkers.txt"))
    adata.uns['operations']=np.append(adata.uns['operations'],inspect.stack()[0][3])
    if save:
        save_adata(adata,sc.settings.figdir)
    return(adata)

#Calculate latent dirichlet allocation of a dataset
def sc_lda(adata,n_components=12,topic_word_prior=None,doc_topic_prior=None,save=False):
    import sklearn
    from sklearn import decomposition
    def moving_average(data,window_width):
        cumsum_vec = np.cumsum(np.insert(data, 0, 0))
        ma_vec = (cumsum_vec[window_width:] - cumsum_vec[:-window_width]) / window_width
        return(ma_vec)
    def table_top_words(word_topic, feature_names):
        return(pd.DataFrame(np.stack([feature_names[(-word_topic).argsort(axis=0)],word_topic[(-word_topic).argsort(axis=0)][:,:,0]],2).reshape(word_topic.shape[0],-1)))

    adata=adata[adata.obs.n_counts.argsort()[::-1],:]
    ldaM=sklearn.decomposition.LatentDirichletAllocation(n_components=n_components,topic_word_prior=topic_word_prior,doc_topic_prior=doc_topic_prior,learning_method="online",verbose=1,n_jobs=-1)
    ldaM.fit(adata.X)
    doc_topic=ldaM.transform(adata.X)
    adata.uns['lda_model']=ldaM
    for i in range(doc_topic.shape[1]):
        adata.obs['lda_'+str(i)]=doc_topic[:,i]

    for i in range(ldaM.components_.shape[0]):
        adata.var['lda_'+str(i)]=np.log10(ldaM.components_[i,:])
    import seaborn
    hm=seaborn.heatmap(doc_topic).get_figure()
    hm.savefig(os.path.join(sc.settings.figdir,"TopicHeatmap.png"))
    plt.clf()
    plt.plot(range(doc_topic.shape[0]),np.std( doc_topic,axis=1))
    plt.title("SD of topic allocation in cells")
    plt.xlabel("RankedCells")
    plt.ylabel("SD")
    plt.savefig(os.path.join(sc.settings.figdir,"TopicCellSD"))
    plt.clf()
    convolvedSD=moving_average(np.std( doc_topic,axis=1),300)
    plt.plot(range(len(convolvedSD)),convolvedSD)
    plt.title("Moving Average Sparsity (SD)")
    plt.xlabel("RankedCells")
    plt.ylabel("SD")
    plt.savefig(os.path.join(sc.settings.figdir,"TopicCellMASD.png"))
    plt.clf()
    plt.plot(range(adata.obs.shape[0]),np.log10(adata.obs.n_counts))
    plt.title("Log Counts per Cell")
    plt.xlabel("RankedCells")
    plt.title("Log Reads")
    plt.savefig(os.path.join(sc.settings.figdir,"LogCounts.png"))
    plt.clf()
    for i in range(doc_topic.shape[1]):
        convolved=moving_average(doc_topic[:,i],300)
        plt.plot(range(len(convolved)),convolved)
    plt.xlabel("RankedCells")
    plt.ylabel("Moving Avg Topic Probability")
    plt.savefig(os.path.join(sc.settings.figdir,"TopicMovingAvg.png"))
    plt.clf()
    table_top_words(ldaM.components_.T,adata.var.index).to_csv(os.path.join(sc.settings.figdir,"TopicMarkers.txt"))
    adata.uns['operations']=np.append(adata.uns['operations'],inspect.stack()[0][3])
    if save:
        save_adata(adata,sc.settings.figdir)
    return(adata)

#TODO
def rna_velocity(adata,loomfile,basis="tsne",save=True):
    sc.tl.addCleanObsNames(adata)
    adata,vdata=sc.tl.rna_velocity(adata,loomfile=loomfile,basis=basis)
    sc.tl.plot_velocity_arrows(adata,basis=basis,cluster='louvain',cluster_colors='louvain_colors')
    plt.savefig(os.path.expanduser(sc.settings.figdir +basis+"_VelocityArrows.pdf"))
    adata.uns['operations']=np.append(adata.uns['operations'],inspect.stack()[0][3])
    if save:
        save_adata(adata,outpath)
        vdata.write(os.path.join(outpath,"_".join([adata.uns['name']]+list(adata.uns['operations']))+"_Velocitydu.h5ad"))
    return(adata,vdata)

#From UMI tools package
def getKneeEstimate(cell_barcode_counts,
                    expect_cells=False,
                    cell_number=False,
                    plotfile_prefix=None,
                    force=True):
    ''' estimate the number of "true" cell barcodes
    input:
         cell_barcode_counts = dict(key = barcode, value = count)
         expect_cells (optional) = define the expected number of cells
         cell_number (optional) = define number of cell barcodes to accept
         plotfile_prefix = (optional) prefix for plots
    returns:
         List of true barcodes
    '''
    from umi_tools import umi_methods
    from collections import Counter
    from functools import partial
    from scipy.signal import argrelextrema
    from scipy.stats import gaussian_kde
    import umi_tools.Utilities as U
    import matplotlib.lines as mlines
    # very low abundance cell barcodes are filtered out (< 0.001 *
    # the most abundant)
    threshold = 0.001 * cell_barcode_counts.most_common(1)[0][1]

    counts = sorted(cell_barcode_counts.values(), reverse=True)
    counts_thresh = [x for x in counts if x > threshold]
    log_counts = np.log10(counts_thresh)

    # guassian density with hardcoded bw
    density = gaussian_kde(log_counts, bw_method=0.1)

    xx_values = 10000  # how many x values for density plot
    xx = np.linspace(log_counts.min(), log_counts.max(), xx_values)

    local_min = None

    if cell_number:  # we have a prior hard expectation on the number of cells
        threshold = counts[cell_number]

    else:
        local_mins = argrelextrema(density(xx), np.less)[0]
        local_mins_counts = []

        for poss_local_min in local_mins[::-1]:

            passing_threshold = sum([y > np.power(10, xx[poss_local_min])
                                     for x, y in cell_barcode_counts.items()])
            local_mins_counts.append(passing_threshold)

            if not local_min:   # if we have selected a local min yet
                if expect_cells:  # we have a "soft" expectation
                    if (passing_threshold > expect_cells * 0.01 and
                        passing_threshold <= expect_cells):
                        local_min = poss_local_min

                else:  # we have no prior expectation
                    # TS: In abscence of any expectation (either hard or soft),
                    # this set of heuristic thresholds are used to decide
                    # which local minimum to select.
                    # This is very unlikely to be the best way to achieve this!
                    if (poss_local_min >= 0.2 * xx_values and
                        (log_counts.max() - xx[poss_local_min] > 0.5 or
                         xx[poss_local_min] < log_counts.max()/2)):
                        local_min = poss_local_min

        if local_min is None and force:
            local_min = min(local_mins)

        if local_min is not None:
            threshold = np.power(10, xx[local_min])

    if cell_number or local_min is not None:
        final_barcodes = set([
            x for x, y in cell_barcode_counts.items() if y > threshold])
    else:
        final_barcodes = None

    if plotfile_prefix:

        # colour-blind friendly colours - https://gist.github.com/thriveth/8560036
        CB_color_cycle = ['#377eb8', '#ff7f00', '#4daf4a',
                          '#f781bf', '#a65628', '#984ea3',
                          '#999999', '#e41a1c', '#dede00']
        user_line = mlines.Line2D(
            [], [], color=CB_color_cycle[0], ls="dashed",
            markersize=15, label='User-defined')
        selected_line = mlines.Line2D(
            [], [], color=CB_color_cycle[0], ls="dashed", markersize=15, label='Selected')
        rejected_line = mlines.Line2D(
            [], [], color=CB_color_cycle[3], ls="dashed", markersize=15, label='Rejected')

        # make density plot
        fig = plt.figure()
        fig1 = fig.add_subplot(111)
        fig1.plot(xx, density(xx), 'k')
        fig1.set_xlabel("Count per cell (log10)")
        fig1.set_ylabel("Density")

        if cell_number:
            fig1.axvline(np.log10(threshold), ls="dashed", color=CB_color_cycle[0])
            lgd = fig1.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.,
                              handles=[user_line],
                              title="Cell threshold")

        elif local_min is None:  # no local_min was accepted
            for pos in xx[local_mins]:
                fig1.axvline(x=pos, ls="dashed", color=CB_color_cycle[3])
            lgd = fig1.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.,
                              handles=[selected_line, rejected_line],
                              title="Possible thresholds")
        else:
            for pos in xx[local_mins]:
                if pos == xx[local_min]:  # selected local minima
                    fig1.axvline(x=xx[local_min], ls="dashed", color=CB_color_cycle[0])
                else:
                    fig1.axvline(x=pos, ls="dashed", color=CB_color_cycle[3])

            lgd = fig1.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.,
                              handles=[selected_line, rejected_line],
                              title="Possible thresholds")

        fig.savefig("%s_cell_barcode_count_density.png" % plotfile_prefix,
                    bbox_extra_artists=(lgd,), bbox_inches='tight')

        # make knee plot
        fig = plt.figure()
        fig2 = fig.add_subplot(111)
        fig2.plot(range(0, len(counts)), np.cumsum(counts), c="black")

        xmax = len(counts)
        if local_min is not None:
            # reasonable maximum x-axis value
            xmax = min(len(final_barcodes) * 5, xmax)

        fig2.set_xlim((0 - (0.01 * xmax), xmax))
        fig2.set_xlabel("Rank")
        fig2.set_ylabel("Cumulative count")

        if cell_number:
            fig2.axvline(x=cell_number, ls="dashed", color=CB_color_cycle[0])
            lgd = fig2.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.,
                              handles=[user_line],
                              title="Cell threshold")

        elif local_min is None:  # no local_min was accepted
            for local_mins_count in local_mins_counts:
                fig2.axvline(x=local_mins_count, ls="dashed",
                             color=CB_color_cycle[3])
            lgd = fig2.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.,
                              handles=[selected_line, rejected_line],
                              title="Possible thresholds")

        else:
            for local_mins_count in local_mins_counts:
                if local_mins_count == len(final_barcodes):  # selected local minima
                    fig2.axvline(x=local_mins_count, ls="dashed",
                                 color=CB_color_cycle[0])
                else:
                    fig2.axvline(x=local_mins_count, ls="dashed",
                                 color=CB_color_cycle[3])

            lgd = fig2.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.,
                              handles=[selected_line, rejected_line],
                              title="Possible thresholds")

        fig.savefig("%s_cell_barcode_knee.png" % plotfile_prefix,
                    bbox_extra_artists=(lgd,), bbox_inches='tight')

        if local_min is not None:
            colours_selected = [CB_color_cycle[0] for x in range(0, len(final_barcodes))]
            colours_rejected = ["black" for x in range(0, len(counts)-len(final_barcodes))]
            colours = colours_selected + colours_rejected
        else:
            colours = ["black" for x in range(0, len(counts))]

        fig = plt.figure()
        fig3 = fig.add_subplot(111)
        fig3.scatter(x=range(1, len(counts)+1), y=counts,
                     c=colours, s=10, linewidths=0)
        fig3.loglog()
        fig3.set_xlim(0, len(counts)*1.25)
        fig3.set_xlabel('Barcode index')
        fig3.set_ylabel('Count')

        if cell_number:
            fig3.axvline(x=cell_number, ls="dashed", color=CB_color_cycle[0])
            lgd = fig3.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.,
                              handles=[user_line],
                              title="Cell threshold")
        elif local_min is None:  # no local_min was accepted
            for local_mins_count in local_mins_counts:
                fig3.axvline(x=local_mins_count, ls="dashed",
                             color=CB_color_cycle[3])
            lgd = fig3.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.,
                              handles=[selected_line, rejected_line],
                              title="Possible thresholds")
        else:
            for local_mins_count in local_mins_counts:
                if local_mins_count == len(final_barcodes):  # selected local minima
                    fig3.axvline(x=local_mins_count, ls="dashed",
                                 color=CB_color_cycle[0])
                else:
                    fig3.axvline(x=local_mins_count, ls="dashed",
                                 color=CB_color_cycle[3])

            lgd = fig3.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.,
                              handles=[selected_line, rejected_line],
                              title="Possible thresholds")

        fig.savefig("%s_cell_barcode_counts.png" % plotfile_prefix,
                    bbox_extra_artists=(lgd,), bbox_inches='tight')

        if not cell_number:
            with U.openFile("%s_cell_thresholds.tsv" % plotfile_prefix, "w") as outf:
                outf.write("count\taction\n")
                for local_mins_count in local_mins_counts:
                    if local_min and local_mins_count == len(final_barcodes):
                        threshold_type = "Selected"
                    else:
                        threshold_type = "Rejected"

                    outf.write("%s\t%s\n" % (local_mins_count, threshold_type))

    return final_barcodes
