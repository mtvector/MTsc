import scvelo as scv
scv.settings.set_figure_params('scvelo')
import scanpy.api as sc
sc.settings.autoshow=False
sc.settings.autosave=True
sc.settings.figdir='/scrapp2/mtschmitz/data/Exonic/fig'
adata = sc.read_10x_mtx('/scrapp2/mtschmitz/data/Exonic/E40_motor_Out/outs/filtered_gene_bc_matrices/refdata-celranger-mmul8-toplevel/', cache=True)
ldata = scv.read('/scrapp2/mtschmitz/data/Exonic/E40_motor_Out_velocyto/possorted_genome_bam_RWRQ2.loom', cache=True)
adata.var_names_make_unique()
ldata.var_names_make_unique()
adata = scv.utils.merge(adata, ldata)
adata.var_names_make_unique()
print('norm')
scv.pp.filter_genes(adata)
scv.pp.normalize_per_cell(adata)
scv.pp.filter_genes_dispersion(adata)
scv.pp.log1p(adata)
print(adata)
print('moment')
scv.pp.moments(adata, n_pcs=30, n_neighbors=30)
print('velo')
scv.tl.umap(adata)
scv.tl.velocity(adata)
print('graph')
scv.tl.velocity_graph(adata)
scv.tl.velocity_embedding(adata, basis='umap')
scv.pl.velocity_embedding(adata, basis='umap',save='Embed')
scv.pl.velocity_embedding_grid(adata, basis='umap',save='Grid')
scv.pl.velocity_embedding_stream(adata, basis='umap',save='stream')
sc.tl.leiden(adata)

def cell_cycle_score(adata,save=False):
    s_genes=['MCM5', 'PCNA', 'TYMS', 'FEN1', 'MCM2', 'MCM4', 'RRM1', 'UNG', 'GINS2', 'MCM6', 'CDCA7', 'DTL', 'PRIM1', 'UHRF1', 'MLF1IP', 'HELLS', 'RFC2', 'RPA2', 'NASP', 'RAD51AP1', 'GMNN', 'WDR76', 'SLBP', 'CCNE2', 'UBR7', 'POLD3', 'MSH2', 'ATAD2', 'RAD51', 'RRM2', 'CDC45', 'CDC6', 'EXO1', 'TIPIN', 'DSCC1', 'BLM', 'CASP8AP2', 'USP1', 'CLSPN', 'POLA1', 'CHAF1B', 'BRIP1', 'E2F8']
    g2m_genes=['HMGB2', 'CDK1', 'NUSAP1', 'UBE2C', 'BIRC5', 'TPX2', 'TOP2A', 'NDC80', 'CKS2', 'NUF2', 'CKS1B', 'MKI67', 'TMPO', 'CENPF', 'TACC3', 'FAM64A', 'SMC4', 'CCNB2', 'CKAP2L', 'CKAP2', 'AURKB', 'BUB1', 'KIF11', 'ANP32E', 'TUBB4B', 'GTSE1', 'KIF20B', 'HJURP', 'CDCA3', 'HN1', 'CDC20', 'TTK', 'CDC25C', 'KIF2C', 'RANGAP1', 'NCAPD2', 'DLGAP5', 'CDCA2', 'CDCA8', 'ECT2', 'KIF23', 'HMMR', 'AURKA', 'PSRC1', 'ANLN', 'LBR', 'CKAP5', 'CENPE', 'CTCF', 'NEK2', 'G2E3', 'GAS2L3', 'CBX5', 'CENPA']
    sc.tl.score_genes_cell_cycle(adata,s_genes=s_genes,g2m_genes=g2m_genes)
    sc.pl.violin(adata, ['G2M_score','S_score'], groupby='leiden')
    if 'X_tsne' in adata.obsm.keys(): 
        sc.pl.tsne(adata, color=['G2M_score','S_score','phase','leiden'],save='_cc')
    sc.pl.umap(adata, color=['G2M_score','S_score','phase','leiden'],save='_cc')
    return(adata)
cell_cycle_score(adata)
sc.pl.umap(adata, color=['leiden','phase','RBFOX3','SOX2'],save='phases')




