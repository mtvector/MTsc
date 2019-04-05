
tab=read.table("~/Downloads/chromium-shared-sample-indexes-plate.csv",sep = ',',header = F,row.names = 1,stringsAsFactors = F)
rownames(tab)=gsub('SI-GA-','',rownames(tab))
sampletab = read.table('~/Downloads/macaque 10X metadata - Regions (scRNA).tsv',sep = '\t',header=T)
macnames= paste(sampletab[,1],sampletab[,3],sep = '_')
macnames=gsub(' \\(','_',macnames)
macnames=gsub('\\)','',macnames)
macnames=gsub('\\+','_AND_',macnames)
macnames=gsub('10\\%','',macnames)
macnames=gsub(' ','_',macnames)
macnames=gsub('__','_',macnames)
macnames=gsub('-all_plus_choroid','',macnames)

head(transtab)
Reduce(c,t(transtab[1:64,]))[1:10]
#the order
transtab= tab[paste0(sampletab[,'BarcodeRow'],sampletab[,'BarcodeCol']),]
write.table(cbind(rep(macnames,each=4),Reduce(c,t(transtab[1:64,]))),file = '~/Downloads/MacaqueSampleList.txt',sep='\t',quote = F,row.names = F,col.names = F)
