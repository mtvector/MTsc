
macpath='~/Downloads/macaque2/'
metricstab <-  sapply(dir(macpath),function(i){
  print(i)
  read.table(paste0(macpath,'/',i,'/outs/metrics_summary.csv'),header = T,sep = ',',stringsAsFactors = F)
})
rn=rownames(metricstab)
metricstab=apply(metricstab,2,function(x)as.numeric(gsub("%|,","",x)))
rownames(metricstab)=rn

for(i in rownames(metricstab)){
  barplot(metricstab[i,],main=i,las=2,cex.names = .6)
}

plot(metricstab['Number.of.Reads',],metricstab['Mean.Reads.per.Cell',])

metricstab['Mean.Reads.per.Cell',]/metricstab['Number.of.Reads',]
