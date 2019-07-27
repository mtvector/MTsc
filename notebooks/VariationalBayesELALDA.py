#!/usr/bin/env python
# coding: utf-8

# # Ambien't

# In[1]:



import sys, os
import theano
import scipy
from collections import OrderedDict
from copy import deepcopy
import numpy as np
import time
import matplotlib.pyplot as plt
import seaborn as sns
from theano import shared
import theano.tensor as tt
from theano.sandbox.rng_mrg import MRG_RandomStreams
import pymc3 as pm
from pymc3 import math as pmmath
from pymc3 import Dirichlet
from pymc3.distributions import Interpolated
from pymc3.distributions.transforms import t_stick_breaking
plt.style.use('seaborn-darkgrid')
from tqdm import tqdm
from collections import Counter
import scanpy as sc
import pandas as pd
sc.settings.autosave=True
sc.settings.autoshow=False
import anndata
import re

def from_epdf(param, ambient_counts):
    from scipy.stats import gaussian_kde
    smin, smax = np.min(ambient_counts), np.max(ambient_counts)
    width = smax - smin
    x = np.linspace(smin, smax, 200)
    y = gaussian_kde(ambient_counts)(x)
    x = np.concatenate([[-1e9,-100,0], x, [x[-1] + np.log10(10)]])
    y = np.concatenate([[1e-1000,1e-200,1e-2], y, [1e-20]])
    sns.scatterplot(x,y)
    plt.show()
    #return pm.Normal('x', mu=np.mean(ambient_counts), sigma=np.std(ambient_counts))
    return Interpolated(param, x, y)


maxgenes=12000
n_iterations=200
#Number of topics, not including ambient topic
K=10
#adatapath=os.path.expanduser('/scrapp2/mtschmitz/5k_pbmc_v3')
#genomename="refdata-celranger-mmul8-toplevel"
#genomename=''
#samplename='5kPBMC'

headpath='/scrapp2/mtschmitz/macaqueseq2/'
adatapaths=[os.path.join(headpath,x) for x in  os.listdir('/scrapp2/mtschmitz/macaqueseq2/')]
samplenames=[re.sub('_Out','',x) for x in  os.listdir('/scrapp2/mtschmitz/macaqueseq2/')]
#adatapaths=['/scrapp2/mtschmitz/E65-2019A_AND_E65-2019B_MULTI-SEQ_1_Out','/scrapp2/mtschmitz/E65-2019A_AND_E65-2019B_MULTI-SEQ_2_Out','/scrapp2/mtschmitz/E65-2019A_AND_E65-2019B_MULTI-SEQ_3_Out','/scrapp2/mtschmitz/E80-2019_MULTI-SEQ_Out','/scrapp2/mtschmitz/E90-2019_MULTI-SEQ_Out','/scrapp2/mtschmitz/1k_hgmm','/scrapp2/mtschmitz/5k_pbmc_v3','/scrapp2/mtschmitz/E40_motor_Out','E50_motor_Out','E65_motor_Out','E80motor_Out','E100motor_Out']
#genomenames=['','','','','','','',"refdata-celranger-mmul8-toplevel","refdata-celranger-mmul8-toplevel","refdata-celranger-mmul8-toplevel","refdata-celranger-mmul8-toplevel","refdata-celranger-mmul8-toplevel"]
#samplenames=['MS65A','MS65B','MS65C','MS80','MS90','1kHGMM','5kPCMB','E40','E50','E65','E80','E100']
ind=np.random.choice(list(range(len(adatapaths))),size=len(adatapaths)-1)
#ind=[0,1]
for adatapath,samplename in zip(np.array(adatapaths)[ind],np.array(samplenames)[ind]):
    print(samplename, flush=True)
    sc.settings.figdir=os.path.expanduser('~/figs/'+samplename+'')
    if not os.path.exists(sc.settings.figdir):
            os.makedirs(sc.settings.figdir)
    
    if os.path.exists(os.path.join(adatapath,'outs/ambientsubtracted.h5ad')):
        continue
        
    adata = sc.read_10x_mtx(os.path.join(adatapath,'outs/raw_feature_bc_matrix'),cache=True)
    bcs=list(pd.read_csv(os.path.join(adatapath,'outs/filtered_feature_bc_matrix','barcodes.tsv.gz')).iloc[:,0])
    #adata = sc.read_10x_h5('/home/mt/code/data/AlignedOrangutanOrganoid/Exonic/orangutanorganoid_Out/outs/filtered_gene_bc_matrices_h5.h5','refdata-celranger-Pabe2-toplevel')
    sc.pp.filter_genes(adata, min_cells=10,inplace=True)
    sc.pp.filter_cells(adata,min_counts=1,inplace=True)
    plt.close()
    sns.distplot(np.log10(adata.X.sum(1).A1),kde=False)
    plt.savefig(os.path.join(sc.settings.figdir,'DropCountsHist.png'))
    plt.close()
    sc.pp.filter_cells(adata,min_counts=15,inplace=True)
    adata=adata[adata.obs.n_counts.argsort(),:]

    #adata._inplace_subset_obs(np.random.choice(adata.obs.index,10000,replace=False))
    #adata._inplace_subset_var(np.random.choice(adata.var.index,4000,replace=False))
    #sc.pp.filter_genes(adata, min_cells=50,inplace=True)
    #sc.pp.filter_cells(adata,min_counts=5,inplace=True)

    adata.var_names_make_unique()
    freshadata=adata.copy()
    print('original shape',freshadata.shape)
    cell_inds=np.where([x in bcs for x in adata.obs.index])[0]
    junk_inds=np.where([x not in bcs for x in adata.obs.index])[0]

    ambient_counts=np.log10(adata[[x not in bcs for x in adata.obs.index] ,:].X.sum(1).A1)
    adata._inplace_subset_var((-adata[junk_inds,:].X).sum(0).A1.argsort()[0:min(maxgenes,adata.shape[1])])
    phiAmbient = adata[junk_inds,:].X.sum(axis=0)/adata[junk_inds,:].X.sum(axis=1).sum()
    phiAmbientDict=dict(zip(list(adata.var.index), phiAmbient.A1))
    adata._inplace_subset_obs(cell_inds)
    freshadata._inplace_subset_obs(cell_inds)
    
    feature_names=list(adata.var.index)
    tf=adata.X
    n_tokens = np.sum(tf[tf.nonzero()])
    print('Number of tokens in training set = {}'.format(n_tokens))
    print('Sparsity = {}'.format(
        len(tf.nonzero()[0]) / float(tf.shape[0] * tf.shape[1])))

    model1 = pm.Model()
    (D,V)=tf.shape
    print(D,V)
    alpha = np.ones((1, K))
    phi = np.ones((1, V))
    sparse_array=shared(np.array([tf.nonzero()[0],tf.nonzero()[1],tf.data]).T.astype('int32'))
    tt.cast(sparse_array,'int32')
    rowsums=shared(np.sum(tf,axis=1).T)
    sumall=shared(np.sum(tf))
    print("Ambient Mean",np.mean(ambient_counts), flush=True)
    print("Ambient STD",np.std(ambient_counts), flush=True)
    def lognormpdf(mean, sd):
        import math
        def internalnorm(x):
            var = float(sd)**2
            denom = tt.log(2*math.pi*var)*.5
            num = -(x-float(mean))**2/(2*var)
            return(num-denom)
        return internalnorm
    ambient=lognormpdf(np.mean(ambient_counts),np.std(ambient_counts))

    def log_lda(theta, phi,value,rowsums,sumall,phiAmbient=None):
        if phiAmbient is not None:
            phi=tt.concatenate([phi,phiAmbient],axis=0)
        else:
            phi=phi
        ll = value[:,2] * pm.math.logsumexp(tt.log(theta[value[:,0].astype('int32')]+1e-9)+ tt.log(phi.T[value[:,1].astype('int32')]+1e-9),axis=1).ravel()                                                                  
        #ambientll=ambient.distribution.logp(tt.log10((rowsums+1e-9)*(theta[:,theta.shape[1]-1]+1e-9)))
        ambientll=ambient(tt.log10((rowsums+1e-9)*(theta[:,theta.shape[1]-1]+1e-9)))
        #tt.printing.Print('l')(ambient.distribution.logp(tt.log10((rowsums+1e-9)*(theta[:,theta.shape[1]-1]+1e-9))))
        return tt.sum(ll) + tt.sum(ambientll)*(sumall/rowsums.shape[1])
        #return tt.sum(ll)/sumall + tt.sum(ambientll)/rowsums.shape[1]
    
    with model1: 
        theta = pm.Dirichlet("theta", a=alpha, shape=(D, K), transform=t_stick_breaking(1e-9)).astype('float32')
        phi = pm.Dirichlet("phi", a=phi, shape=(K-1, V), transform=t_stick_breaking(1e-9)).astype('float32')
        #ambient=from_epdf('ambient',ambient_counts)
        #doc = pm.DensityDist('doc', log_lda, observed=dict(theta=theta, phi=phi,ambient=ambient, value=sparse_array,phiAmbient=np.matrix([phiAmbientDict[x] for x in feature_names]),rowsums=rowsums,sumall=sumall))
        doc = pm.DensityDist('doc', log_lda, observed=dict(theta=theta,
            phi=phi,
            value=sparse_array,phiAmbient=np.matrix([phiAmbientDict[x] for x in feature_names]),rowsums=rowsums,sumall=sumall))

    eta = .5
    s = shared(eta)
    def reduce_rate(a, h, i):
        s.set_value(eta/((i/2)+1)**.5)    
    with model1:    
        inference = pm.ADVI()
        approx = pm.fit(n=n_iterations,method= inference,obj_optimizer=pm.adam(learning_rate=s),callbacks=[reduce_rate,pm.callbacks.CheckParametersConvergence(diff='absolute')])


    tr1 = approx.sample(draws=1000)
    advi_elbo = pd.DataFrame(
        {'log-ELBO': -np.log(approx.hist),
         'n': np.arange(approx.hist.shape[0])})
    plt.clf()
    sns.lineplot(y='log-ELBO', x='n', data=advi_elbo)
    plt.savefig(os.path.join(sc.settings.figdir,'ELBO.png'))
    theta=tr1['theta'].mean(0)
    theta=anndata.AnnData(theta,var=pd.DataFrame(index=['lda_'+str(i) for i in range(K) ]),obs=pd.DataFrame(index=list(adata.obs.index)))

    # In[15]:
    phi=tr1['phi'].mean(0)
    phi=anndata.AnnData(phi,var=pd.DataFrame(index=feature_names),obs=pd.DataFrame(index=['lda_'+str(i) for i in range(K-1) ]))
    
    #phizeros=set(freshadata.var.index)-set(phiAmbientDict.keys())
    #phiAmbientDict.update(dict.fromkeys(list(phizeros),0))
    last=set(freshadata.var.index)-set(phiAmbientDict.keys())
    for i in range(phi.shape[0]):
        freshadata.var['lda_'+str(i)]=0
        freshadata.var.loc[list(phiAmbientDict.keys()),'lda_'+str(i)]=phi[i,:][:,list(phiAmbientDict.keys())].X
    freshadata.var['lda_'+str(K-1)]=0
    #Resort indices so the ones we're subtracting from are listed first
    freshadata.var.loc[list(phiAmbientDict.keys()),'lda_'+str(K-1)]=list(phiAmbientDict.values())
    fullinds=list(freshadata.var.index)
    [fullinds.insert(0, fullinds.pop(fullinds.index(i))) for i in feature_names[::-1]]
    freshadata=freshadata[:,fullinds]
    #Reassign phiAmbient now that matrix is reordered
    phiAmbientDict=dict(zip(list(freshadata.var.index),list(freshadata.var['lda_'+str(K-1)])))
    #phiAmbient=np.array(freshadata.var['lda_'+str(K-1)])    
    
    tmpmat=scipy.sparse.csr_matrix(freshadata.shape)
    
    for ci,c in tqdm(enumerate(freshadata.obs.index)):
        farow=freshadata[c,:].X
        rowinds=farow.nonzero()[0]
        countz=farow[rowinds]
        rownames=freshadata.var.index[rowinds]
        indmultiset=[[rowinds[i]]*int(countz[i]) for i in range(len(countz))]
        vals = [item for sublist in indmultiset for item in sublist]
        #phimultiset=[[phiAmbientDict[rowname]]*int(countz[i]) for i,rowname in enumerate(rownames)]
        phimultiset=[[phiAmbientDict[rowname]/countz[i]]*int(countz[i]) for i,rowname in enumerate(rownames)]
        pvals = [item for sublist in phimultiset for item in sublist]
        '''vals=[]
        pvals=[]
        for ii,i in zip(freshadata[c,:].X.nonzero()[0],freshadata.var.index[freshadata[c,:].X.nonzero()[0]]):
            #make list of all genes observed in cell, repeated by number of counts
            t0=time.process_time()
            vals=vals+[ii]*int(freshadata[c,i].X)
            print("0",time.process_time()-t0)
            t0=time.process_time()
            pvals=pvals+[phiAmbientDict[i]]*int(freshadata[c,i].X)
            print("1",time.process_time()-t0)'''
        if len(vals)>0:
            #Select counts based on phiAmbient
            cellcount=np.sum(freshadata[c,:][:,feature_names].X)
            countremove=Counter(np.random.choice(vals,replace=False,size=min(len(vals),int(theta[c,theta.shape[1]-1].X*cellcount)),p=np.array(pvals)/sum(pvals))) 
            tmpmat[ci,list(countremove.keys())]=list(countremove.values())

    print(tmpmat.sum(), flush=True)
    #subtract the values that are probabalistically likely to be ambient
    freshadata.X=freshadata.X-tmpmat
    print(np.corrcoef(tmpmat.sum(0),[phiAmbientDict[i] for i in freshadata.var.index]), flush=True)
    
    for i in range(theta.shape[1]):
        freshadata.obs['lda_'+str(i)]=theta[:,i].X
    freshadata._inplace_subset_obs([x<.5 for x in freshadata.obs['lda_'+str(K-1)]])
    print(freshadata.obs.index)
    freshadata.write(os.path.join(adatapath,'outs/ambientsubtracted.h5ad'))
    print('saved subtracted mat')
    
    #Now plot Umaps for stuff
    sc.pp.normalize_per_cell(adata, counts_per_cell_after=1e4)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, n_top_genes=5000)
    sc.pp.scale(adata, max_value=10)
    sc.pp.pca(adata)
    sc.pp.neighbors(adata)
    sc.tl.umap(adata)
    sc.tl.leiden(adata)
    sc.pl.umap(adata, color=['leiden'],save="BeforeLeiden")

    #sc.pl.umap(adata, color=['lda_0','lda_1','lda_2','lda_3','lda_4','lda_5','lda_6','lda_7','lda_8','lda_9','n_counts'],save="BeforeLDA")

    sc.tl.rank_genes_groups(adata, 'leiden', method='logreg')
    result = adata.uns['rank_genes_groups']
    groups = result['names'].dtype.names
    df=pd.DataFrame({group + '_' + key[:1]: result[key][group] for group in groups for key in ['names', 'scores']})
    df.to_csv(os.path.join(sc.settings.figdir,samplename+"BeforeleidenLogRegMarkers.csv"))


    #sc.pp.filter_genes(freshadata,min_counts=2,inplace=True)
    sc.pp.filter_cells(freshadata,min_counts=500,inplace=True)
    sc.pp.normalize_total(freshadata, target_sum=1e4)
    sc.pp.log1p(freshadata)
    sc.pp.highly_variable_genes(freshadata,n_top_genes=5000,inplace=True)
    sc.pp.scale(freshadata, max_value=10)
    sc.pp.pca(freshadata)
    sc.pp.neighbors(freshadata)
    sc.tl.umap(freshadata)
    sc.tl.leiden(freshadata)
    sc.pl.umap(freshadata, color=['leiden'],save="AfterLeiden")

    sc.pl.umap(freshadata, color=['lda_0','lda_1','lda_2','lda_3','lda_4','lda_5','lda_6','lda_7','lda_8','lda_9','n_counts'],save="AfterLDA")


    sc.tl.rank_genes_groups(freshadata, 'leiden', method='logreg')
    result = freshadata.uns['rank_genes_groups']
    groups = result['names'].dtype.names
    df=pd.DataFrame({group + '_' + key[:1]: result[key][group] for group in groups for key in ['names', 'scores']})
    df.to_csv(os.path.join(sc.settings.figdir,samplename+"AfterleidenLogRegMarkers.csv"))
