import matplotlib
matplotlib.use('agg') # plotting backend compatible with screen
import scanpy.api as sc
import pandas as pd
import numpy as np
import logging
import os
import velocyto as vcy

logging.basicConfig(level=logging.INFO)

def downsampleVelocyto(dat, numG, numC):
    dat.A=dat.A[0:numG,0:numC]
    dat.U=dat.U[0:numG,0:numC]
    dat.S=dat.S[0:numG,0:numC]
    for k in dat.ca.keys():
        dat.ca[k]=dat.ca[k][0:numC]
    for k in dat.ra.keys():
        dat.ra[k]=dat.ra[k][0:numG]
    vlm.initial_Ucell_size=vlm.initial_Ucell_size[0:numC]
    vlm.initial_cell_size=vlm.initial_cell_size[0:numC]
    return(dat)

logging.info("Start")

sc.settings.verbosity = 2
sc.settings.autosave = True # save figures, do not show them
sc.settings.set_figure_params(dpi=300) # set sufficiently high resolution for saving

inputfile = os.path.expanduser('/ye/yelabstore2/mtschmitz/seq/AlignedOrangutanOrganoid/Exonic/orangutanorganoid_Out/outs')
# Run louvain clustering on true gene expression values
velocityFile = os.path.expanduser('/ye/yelabstore2/mtschmitz/seq/AlignedOrangutanOrganoid/Exonic/orangutanorganoid_Out/velocyto/orangutanorganoid_Out.loom')
vlm = vcy.VelocytoLoom(velocityFile)

logging.info("Ck2")

vlm.normalize("S", size=True, log=True)
vlm.S_norm  # contains log normalized
vlm.filter_cells(bool_array=vlm.initial_Ucell_size > np.percentile(vlm.initial_Ucell_size, 0.5))
logging.info("Ck3")

vlm.set_clusters(vlm.ca["Clusters"])
# Select genes that expressed above threshold for total number of molecules in total and a minimum number of cells containing said molecule.
vlm.score_detection_levels()
vlm.filter_genes(by_detection_levels=True)
# Perform feature selection
vlm.score_cv_vs_mean(5000, max_expr_avg=35)
vlm.filter_genes(by_cv_vs_mean=True)
logging.info("Ck4")
vlm._normalize_S(relative_size=vlm.S.sum(0),
             target_size=vlm.S.sum(0).mean())
vlm._normalize_U(relative_size=vlm.U.sum(0),
             target_size=vlm.U.sum(0).mean())

MasterCells = [vlm.ca['CellID'][ii][-17:-1] for ii in range(len(vlm.ca['CellID']))]
MasterGenes = vlm.ra['Gene'].tolist()
logging.info('Normalization and filtering complete.')
logging.info('Filtered barcodes identified:')
logging.info(str('Sample barcode: ' + MasterCells[0]))
logging.info(str('Number of barcodes: ' + str(np.shape(MasterCells)[0])))
logging.info('Filtered genes identified')
logging.info(str('Sample gene: ' + MasterGenes[0]))
logging.info(str('Number of genes: ' + str(np.shape(MasterGenes)[0])))

logging.info("Ck5")
vlm.perform_PCA()
logging.info("Ck6")
vlm.knn_imputation(n_pca_dims=20, k=500, balanced=True, b_sight=3000, b_maxl=1500, n_jobs=16)
logging.info("Ck7")
vlm.fit_gammas()
vlm.plot_phase_portraits(["SOX2", "TBR1"])
logging.info("Portraits")
matplotlib.savefig(os.path.expanduser('~/OrangutanOrganoidPhasePortraits.pdf'))

vlm.predict_U()
vlm.calculate_velocity()
vlm.calculate_shift(assumption="constant_velocity")
vlm.extrapolate_cell_at_t(delta_t=1.)

from sklearn.manifold import TSNE
bh_tsne = TSNE()
vlm.ts = bh_tsne.fit_transform(vlm.pcs[:, :25])

vlm.estimate_transition_prob(hidim="Sx_sz", embed="ts", transform="sqrt", psc=1,
                             n_neighbors=3500, knn_random=True, sampled_fraction=0.5)
vlm.calculate_embedding_shift(sigma_corr = 0.05, expression_scaling=True)

vlm.calculate_grid_arrows(smooth=0.8, steps=(40, 40), n_neighbors=300)
plt.figure(None,(20,10))
logging.info("TSNE")
vlm.plot_grid_arrows(quiver_scale=0.6,
                    scatter_kwargs_dict={"alpha":0.35, "lw":0.35, "edgecolor":"0.4", "s":38, "rasterized":True}, min_mass=24, angles='xy', scale_units='xy',
                    headaxislength=2.75, headlength=5, headwidth=4.8, minlength=1.5,
                    plot_random=True, scale_type="absolute")
matplotlib.savefig(os.path.expanduser('~/OrangutanOrganoidArrowTSNE.pdf'))
vlm.to_hdf5(os.path.expanduser('/ye/yelabstore2/mtschmitz/seq/AlignedOrangutanOrganoid/Exonic/orangutanorganoid_Out/velocyto/orangutanorganoid_Out_analysis.hdf5'))
