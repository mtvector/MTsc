if (!requireNamespace("BiocManager", quietly=TRUE)){
  install.packages("BiocManager")}
#BiocManager::install(version = "3.9")
#BiocManager::install("Seurat")
#BiocManager::install("celda")
#Arg should be path to 10x cellranger 3.0+ output directory. Should be trailing slash
args = commandArgs(trailingOnly=TRUE)
library(celda)
library(Seurat)
#hgmm=Read10X("/scrapp2/mtschmitz/macaqueseq2/E65-2019A_AND_E65-2019B_MULTI-SEQ_1_Out/outs/filtered_feature_bc_matrix/")
hgmm=Read10X(paste0(args[1],'outs/filtered_feature_bc_matrix/'))

#clusters=read.table('/scrapp2/mtschmitz/macaqueseq2/E65-2019A_AND_E65-2019B_MULTI-SEQ_1_Out/outs/analysis/clustering/graphclust/clusters.csv',sep=',',header=T,stringsAsFactors=F)[,2]
clusters=read.table(paste0(args[1],'outs/analysis/clustering/graphclust/clusters.csv'),header=T,stringsAsFactors=F,sep=',')[,2]
decontxModel= decontX(counts = as.matrix(hgmm),z=clusters)
print(decontxModel$resList$logLikelihood)
print(dim(decontxModel$resList$estNativeCounts))
mat=CreateSeuratObject(decontxModel$resList$estNativeCounts)
print(mat)
FindVariableFeatures(mat,nfeatures = nrow(decontxModel$resList$estNativeCounts))
print(mat)
as.loom(mat, filename = paste0(args[1],"outs/deconted.loom"), verbose = FALSE)
