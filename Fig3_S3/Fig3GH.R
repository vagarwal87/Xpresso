library(ROCR)

a=read.delim("all_crossvalidated_predictions_mESC.txt")
colnames(a)[2:3]=c("mESCPred","mESCActual")
c=read.delim("ensembl2geneName_v90_mm10.txt")
colnames(c)[2]="geneName"
a=merge(a,c,by=1,all.x=T)

a$PredAdj=predict(lm(mESCActual~mESCPred, a))
a$resid = a$mESCActual-a$PredAdj
nrow(a)
a = a[a$mESCActual > min(a$mESCActual)+1,]
nrow(a)
miR = fastread("zcat Summary_Counts.default_predictions.txt.gz")
nrow(a)

values = sapply(unique(miR$"miRNA family"), function(fam){
	miR2 = miR[miR$"Species ID"==10090 & miR$"miRNA family"==fam,c(2,3,16)]
	merged = merge(miR2, a, by.x=1,by.y=4, all.y=T)
	merged[is.na(merged)]=0
	
	if (sum(merged[,"Cumulative weighted context++ score"] != 0) > 10){
		c(0,cor(merged$resid, as.numeric(merged[,"Cumulative weighted context++ score"]), method='spearman'))
	}
	else{
		c(0,0)
	}
})
values=t(values)

pdf("Fig3GH.pdf", width=8, height=8)
par(mar=c(7,7,5,5), mgp = c(5, 1.5, 0))
miRs=read.delim("mouseESC_GSE76288_miRNA_counts_Denzler.txt")
miRs=miRs[,c("miRNA_seed","Embryonic_stem_cells._Average_RPM")]
miRs$miRNA_seed=gsub("T","U",miRs$miRNA_seed)
miRs=aggregate(miRs$Embryonic_stem_cells._Average_RPM,by=list(miRs$miRNA_seed), sum)
miRs=miRs[order(miRs$x, decreasing=T),]
miRs2=rbind(miRs[1:10,], c("Other", sum(miRs[11:nrow(miRs),"x"])))
pie(as.numeric(miRs2[,2]), labels = miRs2[,1],col=c('orange','blue','red','purple','darkolivegreen1','magenta','brown', 'cyan', 'yellow','grey','black'), main=paste("Top 10 miRNA families in mESCs"), clockwise = T, cex.main = 2, cex=1.8)

miRs=merge(values, miRs, by.x=0, by.y=1, all.x=T)
miRs$color="black"
miRs$color[miRs$Row.names %in% miRs2[,1]]="red"
colnames(miRs)[3]="spearman"
p1 <- hist(miRs$spearman[miRs$color!="red"], 50, plot=F)
p2 <- hist(miRs$spearman[miRs$color=="red"], 50, plot=F)
plot( p1, col=rgb(0,0,1,1/2), xlim=c(-0.02,0.06))
plot( p2, col=rgb(1,0,0,1/2), xlim=c(-0.02,0.06), add=T)
dev.off()