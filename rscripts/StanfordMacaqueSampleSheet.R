sampletab = read.table('~/stanfordseqindices.csv',sep = ',',header=T,stringsAsFactors = F)
lanes=rep(c(1,2,3,4),nrow(sampletab))
sampletab=sampletab[,-1]
tinytab=cbind(Lane=lanes,Sample=sampletab[,1],Index=sampletab[,ncol(sampletab)])



name10x=sapply(sampletab[,2],function(s){
  rownames(tab)[which(tab[,1]==s)]
})
tinytab=cbind(Lane=lanes,Sample=sampletab[,1],Index=name10x)
#tinytab=cbind(tinytab,Sample_Project=rep('2018MacaqueStanford',nrow(tinytab)))
write.table(tinytab,file = '~/Downloads/stanford_tiny_indexsheet.txt',quote=F,sep=',',row.names = F)

sampletab=t(rbind(lanes,t(rep(sampletab[,1],each=4)),as.vector(t(sampletab[,2:5]))))

tab=read.table("~/Downloads/chromium-shared-sample-indexes-plate.csv",sep = ',',header = F,row.names = 1,stringsAsFactors = F)


sampletab[,6]="ya"
?read.table

