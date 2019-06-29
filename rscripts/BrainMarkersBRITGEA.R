library(BRETIGEA)
df <- BRETIGEA::markers_df_human_brain
table(df[,2])
df[,2] <- gsub('ast','Astrocyte', df[,2])
df[,2] <- gsub('end','Endothelial Cell', df[,2])
df[,2] <- gsub('mic','Microglia', df[,2])
df[,2] <- gsub('neu','Neuron', df[,2])
df[,2] <- gsub('oli','Oligodendrocyte', df[,2])
df[,2] <- gsub('opc','Oligodendrocyte Progenitor Cell', df[,2])
table(df[,2])

markers <- read.csv2('https://docs.google.com/spreadsheets/d/e/2PACX-1vTz5a6QncpOOO-f3FHW2Edomn7YM5mOJu4z_y07OE3Q4TzcRr14iZuVyXWHv8rQuejzhhPlEBBH1y0V/pub?gid=1154528422&single=true&output=tsv',sep = '\t')
markers

cellTypeSpecific <- lapply(unique(df[,2]),function(ct){
  subdf <- df[df[,2]==ct,]
  subdf <- subdf[!(df[,1]%in%markers[,1]),]
  subdf[1:20,]
})

final <- Reduce(rbind,cellTypeSpecific)
finalFormatted <- cbind(final[,1],'',final[,2],'','BRETIGEA_top20')
write.table(finalFormatted,file = '~/markers/BRETIGEA.txt',sep = "\t",row.names = F,quote = F)
