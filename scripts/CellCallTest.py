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
from sklearn import linear_model
from sklearn import neighbors
from sklearn import semi_supervised
from sklearn import ensemble
from matplotlib import pyplot as plot
from matplotlib.backends.backend_pdf import PdfPages
from umi_tools import umi_methods
from collections import Counter
import logging as logg

logg.basicConfig(level=logg.INFO)

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
    import matplotlib.lines as mlines
    from functools import partial
    from scipy.signal import argrelextrema
    from scipy.stats import gaussian_kde
    import umi_tools.Utilities as U
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
            local_min = max(local_mins)

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

def filterSupervisedPCA(adata,lmLR=sklearn.linear_model.LogisticRegression(),knee=True,expect_cells=None):
    name=adata.name+lmLR.__class__.__name__+str(knee)
    adata=adata.copy()
    which = lambda lst:list(np.where(lst)[0])
    plt.switch_backend('agg')
    if knee:
        knee_cells=list(getKneeEstimate(Counter(adata.obs.n_counts.to_dict()),expect_cells=expect_cells,plotfile_prefix=name))
        inthresh=np.array([x in knee_cells for x in adata.obs.index])
        countThreshHigh=np.min(adata.obs.n_counts[knee_cells])
    else:
        countThreshHigh=adata.obs.n_counts.quantile(.99)*.2
        inthresh=np.array(adata.obs['n_counts'])>countThreshHigh    
    inthresh[np.random.choice(which(np.invert(inthresh)), size= min(np.sum(inthresh),np.sum(np.invert(inthresh))), replace=False)]=True   
    obs_copy=adata.obs
    sc.pp.filter_genes(adata, min_cells=50)
    sc.pp.log1p(adata)
    sc.pp.scale(adata, max_value=10)
    sc.pp.pca(adata)
    obs_copy=pd.concat([ obs_copy.reset_index(drop=True),pd.DataFrame(adata.obsm['X_pca'][:,0:9]).reset_index(drop=True)], axis=1)
    obs_copy_sub=obs_copy.loc[inthresh,:]
    trainTable=obs_copy_sub.iloc[:,1:]
    trainTable.n_counts=np.log10(trainTable.n_counts)
    lmLR.fit(X=trainTable,y=trainTable.n_counts>np.log10(countThreshHigh))
    plt.clf()
    plt.hist(np.log10(obs_copy['n_counts']),bins=100)
    plt.savefig(os.path.join(sc.settings.figdir,name+"_4.png"))
    plt.hist(trainTable['n_counts'],bins=100)
    plt.savefig(os.path.join(sc.settings.figdir,name+"_5.png"))
    obs_copy.n_counts=np.log10(obs_copy.n_counts)
    #adata.obs['p_cell'] =lmLR.predict_proba(obs_copy.iloc[:,1:])[:,1]
    predict=lmLR.predict(obs_copy.iloc[:,1:])
    adata.obs['c_cell'] =list(map(str,predict))
    print(adata.obs.c_cell.value_counts())
    sc.pl.scatter(adata,color='c_cell', x='n_counts', y='percent_mito',title=name,save=name+"_1.pdf")
    sc.pl.scatter(adata,color='c_cell', x='n_counts', y='percent_ribo',title=name,save=name+"_2.pdf")
    sc.pl.scatter(adata,color='c_cell', x='n_counts', y='n_genes',title=name,save=name+"_3.pdf")
    return([name]+list(adata.obs.c_cell.value_counts())+[countThreshHigh])

def adataToTrainState(filename,refname):
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
    scanpy.api.pp.filter_genes_dispersion(adata,n_top_genes=np.sum(np.sum(adata.X, axis=0)>0))
    return(adata)

refname="refdata-celranger-mmul8-toplevel"
    
Qlist=[True,False]
    
sc.settings.autosave=True
sc.settings.autoshow=False

basepath=os.path.expanduser("~/tmp/Macaque/")
fileList = ["E40_MGE_Out","E40_V1_Out","E100temporal_Out_180326_A00269_0056_AH7Y2KDMXX-E100temporal_S16","E100hippo_Out",
           "E40_DRG_Out","E40_Pre-optic_Out","orangutanorganoid_Mmul8_Out"]

fileList=os.listdir(os.path.join(basepath,'Exonic'))

modelList = [sklearn.neighbors.KNeighborsClassifier(n_neighbors=10),
            sklearn.neighbors.KNeighborsClassifier(n_neighbors=100),
            sklearn.linear_model.LogisticRegression(),
            sklearn.semi_supervised.LabelPropagation(kernel='knn',n_neighbors=100),
            sklearn.ensemble.RandomForestClassifier(n_estimators=100)]    

for f in fileList:
    resultList = []
    directory=os.path.join(basepath,'CellCallingFigs',f)
    adata=adataToTrainState(os.path.join(basepath,'Exonic',f,'outs/raw_gene_bc_matrices_h5.h5'),refname)
    sc.settings.figdir=directory
    if not os.path.exists(directory):
        os.makedirs(directory)
    logg.info(os.path.join(basepath,'Exonic',f,'/outs/raw_gene_bc_matrices_h5.h5'))
    for m in modelList:
        #modelOut=pd.DataFrame(columns=['name','False','True','countThreshStart'])
        modelOut=pd.DataFrame()
        for q in Qlist:
            outlist=filterSupervisedPCA(adata,lmLR=m,knee=q,expect_cells=30000)
            logg.info(outlist)
            modelOut=modelOut.append(outlist)
        resultList.append(modelOut)
    pd.concat(resultList).to_csv(basepath+f+'OutTable'+'.txt')
    

