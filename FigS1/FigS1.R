library(gplots)

x=read.delim(gzfile("57epigenomes.RPKM.pc.gz"), row.names=1)
x$E000=NULL
names=read.delim("EG.name.txt",F)
colnames(x)=paste(colnames(x), gsub("_", " ", as.character(unlist(sapply(colnames(x), function(x) names[names$V1==x, "V2"])))))

y=as.matrix(cor(x, method='spearman'))
colnames(y)=colnames(x)

pdf("FigS1.pdf", height=10, width=10)
par(oma=c(16,1,1,14))
heatmap.2(y, trace="none", breaks=seq(0,1,0.05), #,density.info="none"
symkey=FALSE, cexRow=0.6, cexCol=0.6, dendrogram="row", key=TRUE, col=matlab::jet.colors(20), denscol="black")
dev.off()
