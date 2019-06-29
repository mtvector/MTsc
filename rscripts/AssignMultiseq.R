#devtools::install_github('chris-mcginnis-ucsf/MULTI-seq')
#Need to add symlink to /bin/tar from /bin/gtar
library(ggplot2)
library(deMULTIplex)
library(matrixStats)
for( a in dir('~/Downloads/MULTIseqOutputs/')){
  if(!grepl('barTable.txt',a)){
    next
  }
  bar.table.full=read.table(paste0('~/Downloads/MULTIseqOutputs/',a),header=T,row.names=1,sep='\t',stringsAsFactors = F)
  print(a)
  rownames(bar.table.full)=make.unique(rownames(bar.table.full))
  d=gsub('_barTable.txt','',a)
  
  ## Note: Exclude columns 97:98 (assuming 96 barcodes were used) which provide total barcode UMI counts for each cell. 
  
  bar.table <- bar.table.full[,1:96]
  bar.tsne <- barTSNE(bar.table) 
  pdf(paste0("~/Downloads/MULTIseqOutputs/",d,"bc.check.pdf"))
  hist(rowMaxs(as.matrix(bar.table))/rowSums(bar.table))
  for (i in 3:ncol(bar.tsne)) {
    g <- ggplot(bar.tsne, aes(x = TSNE1, y = TSNE2, color = bar.tsne[,i])) +
      geom_point() +
      scale_color_gradient(low = "black", high = "red") +
      ggtitle(colnames(bar.tsne)[i]) +
      theme(legend.position = "none")
    print(g)
  }
  ## Round 1 -----------------------------------------------------------------------------------------------------
  ## Perform Quantile Sweep
  bar.table_sweep.list <- list()
  n <- 0
  for (q in seq(0.01, 0.99, by=0.02)) {
    print(q)
    n <- n + 1
    bar.table_sweep.list[[n]] <- classifyCells(bar.table, q=q)
    names(bar.table_sweep.list)[n] <- paste("q=",q,sep="")
  }
  
  ## Identify ideal inter-maxima quantile to set barcode-specific thresholds
  threshold.results1 <- findThresh(call.list=bar.table_sweep.list)
  ggplot(data=threshold.results1$res, aes(x=q, y=Proportion, color=Subset)) + geom_line() + theme(legend.position = "none") +
    geom_vline(xintercept=threshold.results1$extrema, lty=2) + scale_color_manual(values=c("red","black","blue"))
  
  ## Finalize round 1 classifications, remove negative cells
  round1.calls <- classifyCells(bar.table, q=findQ(threshold.results1$res, threshold.results1$extrema))
  neg.cells <- names(round1.calls)[round1.calls == "Negative"]
  bar.table <- bar.table[!rownames(bar.table) %in% neg.cells, ]
  
  ## Round 2 -----------------------------------------------------------------------------------------------------
  bar.table_sweep.list <- list()
  n <- 0
  for (q in seq(0.01, 0.99, by=0.02)) {
    print(q)
    n <- n + 1
    bar.table_sweep.list[[n]] <- classifyCells(bar.table, q=q)
    names(bar.table_sweep.list)[n] <- paste("q=",q,sep="")
  }
  
  threshold.results2 <- findThresh(call.list=bar.table_sweep.list)
  round2.calls <- classifyCells(bar.table, q=findQ(threshold.results1$res, threshold.results2$extrema))
  neg.cells <- c(neg.cells, names(round2.calls)[round2.calls == "Negative"])
  
  ## Repeat until all no negative cells remain (usually 3 rounds)...
  final.calls <- c(round2.calls[round2.calls != "Negative"], rep("Negative",length(neg.cells)))
  names(final.calls) <- c(names(round2.calls)[round2.calls != "Negative"],neg.cells)
  print(d)
  write.table(sort(table(final.calls),decreasing = T),file=paste0('~/Downloads/MULTIseqOutputs/',d,'_bctable.txt'),sep='\t',quote=F,col.names=T,row.names = T)
  # reclass.cells <- findReclassCells(as.matrix(bar.table.full[,1:96]), names(final.calls)[which(final.calls=="Negative")])
  # reclass.res <- rescueCells(bar.table.full, final.calls, reclass.cells)  
  # 
  # ggplot(reclass.res[-1, ], aes(x=ClassStability, y=MatchRate_mean)) + 
  #   geom_point() + xlim(c(nrow(pool.reclass.res)-1,1)) + 
  #   ylim(c(0,1.05)) +
  #   geom_errorbar(aes(ymin=MatchRate_mean-MatchRate_sd, ymax=MatchRate_mean+MatchRate_sd), width=.1) +
  #   geom_hline(yintercept = reclass.res$MatchRate_mean[1], color="red") +
  #   geom_hline(yintercept = reclass.res$MatchRate_mean[1]+3*reclass.res$MatchRate_sd[1], color="red",lty=2) +
  #   geom_hline(yintercept = reclass.res$MatchRate_mean[1]-3*reclass.res$MatchRate_sd[1], color="red",lty=2)
  dev.off()
  # final.calls.rescued <- final.calls
  # rescue.ind <- which(reclass.cells$ClassStability >= 16) ## Note: Value will be dataset-specific
  # final.calls.rescued[rownames(reclass.cells)[rescue.ind]] <- reclass.cells$Reclassification[rescue.ind]
  write.table(t(as.matrix(final.calls)),file=paste0('~/Downloads/MULTIseqOutputs/',d,'_finalcalls.txt'),sep='\t',quote=F,col.names=T,row.names = T)
}

