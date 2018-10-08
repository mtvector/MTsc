import matplotlib
matplotlib.use('agg') # plotting backend compatible with screen
import scanpy.api as sc
import pandas as pd
import numpy as np
import logging
import os
import velocyto as vcy

logging.basicConfig(level=logging.INFO)


logging.info("Start")

sc.settings.verbosity = 2
sc.settings.autosave = True # save figures, do not show them
sc.settings.set_figure_params(dpi=300) # set sufficiently high resolution for saving

inputfile = os.path.expanduser('/ye/yelabstore2/mtschmitz/seq/AlignedOrangutanOrganoid/Exonic/orangutanorganoid_Out/outs')
# Run louvain clustering on true gene expression values
velocityFile = os.path.expanduser('/ye/yelabstore2/mtschmitz/seq/AlignedOrangutanOrganoid/Exonic/orangutanorganoid_Out/velocyto/orangutanorganoid_Out.loom')

adata = sc.read_10x_h5('/home/mt/code/data/AlignedOrangutanOrganoid/Exonic/orangutanorganoid_Out/outs/filtered_gene_bc_matrices_h5.h5','refdata-celranger-Pabe2-toplevel')
adata=sc.tl.rna_velocity(adata,os.path.expanduser('~/code/data/AlignedOrangutanOrganoid/Exonic/orangutanorganoid_Out/velocyto/orangutanorganoid_Out.loom'))
