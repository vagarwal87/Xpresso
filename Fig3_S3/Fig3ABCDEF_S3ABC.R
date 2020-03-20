#### MOUSE ######

a=read.delim("all_crossvalidated_predictions_mESC.txt")
b=read.delim("all_crossvalidated_predictions_mouse.txt")
colnames(a)[2:3]=c("mESCPred","mESCActual")
a=merge(a,b,by=1)
c=read.delim("ensembl2geneName_v90_mm10.txt")
colnames(c)[2]="geneName"
a=merge(a,c,by=1,all.x=T)
a[a=='']="NA"

summary(lm(mESCActual~mESCPred, a))
summary(lm(mESCActual~Pred, a))

# Table1 Moorthy et al, 2017
a$color='black'
a$color[a$geneName=="Sall1"]='red'
a$color[a$geneName=="Tet1"]='red'
a$color[a$geneName=="Prkcg"]='red'
a$color[a$geneName=="AU018091"]='red'
a$color[a$geneName=="Med13l"]='red'
a$color[a$geneName=="Macf1"]='red'
a$color[a$geneName=="Ranbp17"]='red'
a$color[a$geneName=="Cbfa2t2"]='red'
a$color[a$geneName=="Esrrb"]='red'
a$color[a$geneName=="Dppa5a"]='red'
a$color[a$geneName=="Ooep"]='red'
a$color[a$geneName=="Mcl1"]='red'
a$color[a$geneName=="Etl4"]='red'

# A few genes from Whyte et al, 2013; Dowen et al, 2014; Hnisz et al 2013
a$color[a$geneName=="Pou5f1"]='red' #same as Oct4
a$color[a$geneName=="Sox2"]='red'
a$color[a$geneName=="Nanog"]='red'
a$color[a$geneName=="Klf4"]='red'
a$color[a$geneName=="Tbx3"]='red'
a$color[a$geneName=="Sall4"]='red'
a$color[a$geneName=="Lefty1"]='red'
a$color[a$geneName=="Lefty2"]='red'
a$color[a$geneName=="Utf1"]='red'
a$color[a$geneName=="Phc1"]='red'
a$color[a$geneName=="Nr5a2"]='red'
a$color[a$geneName=="Lrrc2"]='red'
a$color[a$geneName=="Dppa3"]='red'
a$color[a$geneName=="Prdm14"]='red'

a$PredAdj=predict(lm(mESCActual~mESCPred, a))
a$MedianPred=predict(lm(mESCActual~Pred, a))
a$resid = a$mESCActual-a$PredAdj

nrow(a)

pdf("Fig3DEF_S3C.pdf", width=8, height=8)
par(mar=c(7,7,5,5), mgp = c(5, 1.5, 0))

#Fig3D
smoothScatter(a$mESCPred, a$mESCActual, cex.axis=2, cex.lab=2, bty="n", xlab="Predicted Expression Level, mESC model", ylab="mESC expression level (log10)", xlim=c(-1.5, 2), ylim=c(-1, 4), las=1, cex=.5) #, pch=1, col = a$color
abline(0,1, col="red")
text(a[a$color=="red","mESCPred"], a[a$color=="red","mESCActual"], labels = a[a$color=="red","geneName"], offset = 0.5, col="red")
text(1.5, 4, labels = paste("r^2 =", round(cor(a$mESCPred, a$mESCActual)^2,2)), offset = 0.5, col="black")

#FigS3C
plot.ecdf(a$resid[a$color=='black'], xlim=c(-4,4), ylim=c(0,1), verticals= TRUE, do.points = FALSE, col="black")
plot.ecdf(a$resid[a$color=='red'], verticals= TRUE, do.points = FALSE, col="red", add=T)
legend("topleft", bg="white", bty="n", legend = c(paste("non-enhancer-driven genes, n = ", length(a$resid[a$color=='black']), sep=''),
paste("enhancer-driven genes, n = ", length(a$resid[a$color=='red']), sep=''),
paste("P value: ", formatC(ks.test(a$resid[a$color=='red'],a$resid[a$color=='black'],alternative="less")$p.value, digits = 2, format = 'g'), sep='')), text.col = c("black","red", "black"))

#Fig3E
silenced = read.delim("Boyer_et_al_PCG_repressed.txt",F) #from Boyer et al
silenced = merge(silenced, c, by.x=1, by.y=2)
nrow(silenced)
active = read.delim("Whyte_et_al_superenhancers.txt",F) #from Whyte et al
active = merge(active, c, by.x=1, by.y=2)
nrow(active)
a$color='black'
a$color[a$Gene %in% silenced[,2] & !(a$Gene %in% active[,2])]='blue'
a$color[a$Gene %in% active[,2] & !(a$Gene %in% silenced[,2])]='red'
plot.ecdf(a$resid[a$color=='black'], xlim=c(-2,2), ylim=c(0,1), verticals= TRUE, do.points = FALSE, col="black", main="mESC")
plot.ecdf(a$resid[a$color=='blue'], verticals= TRUE, do.points = FALSE, col="blue", add=T)
plot.ecdf(a$resid[a$color=='red'], verticals= TRUE, do.points = FALSE, col="red", add=T)
legend("topleft", bg="white", bty="n", legend = c(paste("Other genes, n = ", length(a$resid[a$color=='black']), sep=''),
paste("PCG-silenced genes, n = ", length(a$resid[a$color=='blue']), sep=''),
paste("Super-enhancer-associated genes, n = ", length(a$resid[a$color=='red']), sep=''),
paste("PCG-silenced vs Black P value: ", formatC(ks.test(a$resid[a$color=='blue'],a$resid[a$color=='black'],alternative="greater")$p.value, digits = 2, format = 'g'), sep=''),  
paste("Super-enhancer vs Black P value: ", formatC(ks.test(a$resid[a$color=='red'],a$resid[a$color=='black'],alternative="less")$p.value, digits = 2, format = 'g'), sep='')), 
text.col = c("black","blue","red","black"))

#Fig3F
halflife = read.delim("Herzog_mESC_half_life.txt") #from Herzog et al
halflife = halflife[,c(4,7)]
colnames(halflife)[2]='half_life'
halflife$half_life=log2(halflife$half_life)
halflife = merge(halflife, c, by.x=1, by.y=2)
a1=merge(halflife, a, by.x=3, by.y=1)
"mESC half lives measured for this many genes:"
nrow(a1)
a1$quintile <- cut(a1$half_life, breaks=quantile(a1$half_life, probs=seq(0,1, by=0.2), na.rm=TRUE), include.lowest=TRUE)
boxplot(a1$resid~a1$quintile,outline=F, cex=1.5, cex.axis=2, cex.lab=2, cex.main=2, las=2, notch=T, col="red")
cor.test(a1$half_life, a1$resid)
cor.test(a1$half_life, a1$resid, method='spearman')

mouse = a

a=read.delim("all_crossvalidated_predictions_mESC.txt",stringsAsFactors=F)
colnames(a)[2:3]=c("mESCPred","mESCActual")
b=read.delim("Ouyang_mESC_RPKM_ensemblID.txt",F)
b$V2=log10(b$V2+0.1)
a=merge(a,b,by=1,all.x=T)
a$mESCActual=a$V2
a$V2=NULL
b=read.delim("all_crossvalidated_predictions_mouse.txt",stringsAsFactors=F)
colnames(b)[2:3]=c("Pred","Actual")
c=read.delim("mouse.median_expr.txt",F)
c$V2=log10(c$V2+0.1)
b=merge(b,c,by=1)
b$Actual=b$V2
b$V2=NULL
a=merge(a,b,by=1,all=T)
c=read.delim("ensembl2geneName_v90_mm10.txt",F,stringsAsFactors=F)
c=c[,c(1,2,4)]
colnames(c)[2:3]=c("geneName","Description")
a=merge(a,c,by=1,all.x=T)
a[a=='']=NA
writefile(cbind(a$Gene, a$geneName, a$Description, round(a$Pred,3), round(a$Actual,3), round(a$mESCPred,3), round(a$mESCActual,3)), "TableS1_mouse.txt", col.names=T)

### HUMAN #####

a=read.delim(gzfile("57epigenomes.RPKM.pc.gz"))
a$E000=NULL
a[,2:ncol(a)]=log10(a[,2:ncol(a)]+0.1)
a$medianExpr=apply(a[,2:ncol(a)], 1, median)

nrow(a)
b=read.delim("all_crossvalidated_predictions.txt")
a=merge(a,b,by=1)
b=read.delim("all_crossvalidated_predictions_K562.txt")
colnames(b)[2:3]=c("K562Pred","K562Actual")
a=merge(a,b,by=1)
b=read.delim("all_crossvalidated_predictions_GM12878.txt")
colnames(b)[2:3]=c("GM12878Pred","GM12878Actual")
a=merge(a,b,by=1)
nrow(a)

c=read.delim("GSE78709_sure23.plasmid.norm.combined.45.55.minus.promoters.bigWigSignal",F) #from van Aresbergen et al
d=read.delim("GSE78709_sure23.plasmid.norm.combined.45.55.plus.promoters.bigWigSignal",F)
c=rbind(c,d)
c$V6=log10(c$V6+0.1)
c=c[,c(1,6)]
colnames(c)[2]='SuRE'
a=merge(a,c,by=1,all.x=T)

c=read.delim("EnsemblID2GeneName.txt",F)
colnames(c)[2:3]=c("geneName","Description")
a=merge(a,c,by=1)
a[a=='']="NA"

pdf("Fig3ABC_S3ABC.pdf", width=8, height=8)
par(mar=c(7,7,5,5), mgp = c(5, 1.5, 0))

a$resid = a$GM12878Actual-a$GM12878Pred

#FigS3B
silenced = read.delim("diHMM/GM12878/H3K27me3_silenced.txt",F) #from Marco et al
nrow(silenced)
active = read.delim("diHMM/GM12878/superenhancer.txt",F)
nrow(active)
a$color='black'
a$color[a$gene_id %in% silenced[,1] & !(a$gene_id %in% active[,1])]='blue'
a$color[a$gene_id %in% active[,1] & !(a$gene_id %in% silenced[,1])]='red'
plot.ecdf(a$resid[a$color=='black'], xlim=c(-2,2), ylim=c(0,1), verticals= TRUE, do.points = FALSE, col="black", main="GM12878")
plot.ecdf(a$resid[a$color=='blue'], verticals= TRUE, do.points = FALSE, col="blue", add=T)
plot.ecdf(a$resid[a$color=='red'], verticals= TRUE, do.points = FALSE, col="red", add=T)
legend("topleft", bg="white", bty="n", legend = c(paste("Other genes, n = ", length(a$resid[a$color=='black']), sep=''),
paste("Silenced genes, n = ", length(a$resid[a$color=='blue']), sep=''),
paste("Stretch-enhancer-associated genes, n = ", length(a$resid[a$color=='red']), sep=''),
paste("Silenced vs Black P value: ", formatC(ks.test(a$resid[a$color=='blue'],a$resid[a$color=='black'],alternative="greater")$p.value, digits = 2, format = 'g'), sep=''),  
paste("Enhancer vs Black P value: ", formatC(ks.test(a$resid[a$color=='red'],a$resid[a$color=='black'],alternative="less")$p.value, digits = 2, format = 'g'), sep='')), 
text.col = c("black","blue","red","black"))

#Fig3B
a$resid = a$K562Actual-a$K562Pred
silenced = read.delim("diHMM/K562/H3K27me3_silenced.txt",F) #from Marco et al
nrow(silenced)
active = read.delim("diHMM/K562/superenhancer.txt",F)
nrow(active)
a$color='black'
a$color[a$gene_id %in% silenced[,1] & !(a$gene_id %in% active[,1])]='blue'
a$color[a$gene_id %in% active[,1] & !(a$gene_id %in% silenced[,1])]='red'
plot.ecdf(a$resid[a$color=='black'], xlim=c(-2,2), ylim=c(0,1), verticals= TRUE, do.points = FALSE, col="black", main="K562")
plot.ecdf(a$resid[a$color=='blue'], verticals= TRUE, do.points = FALSE, col="blue", add=T)
plot.ecdf(a$resid[a$color=='red'], verticals= TRUE, do.points = FALSE, col="red", add=T)
legend("topleft", bg="white", bty="n", legend = c(paste("Other genes, n = ", length(a$resid[a$color=='black']), sep=''),
paste("silenced genes, n = ", length(a$resid[a$color=='blue']), sep=''),
paste("stretch-enhancer-associated genes, n = ", length(a$resid[a$color=='red']), sep=''),
paste("Silenced vs Black P value: ", formatC(ks.test(a$resid[a$color=='blue'],a$resid[a$color=='black'],alternative="greater")$p.value, digits = 2, format = 'g'), sep=''), 
paste("Enhancer vs Black P value: ", formatC(ks.test(a$resid[a$color=='red'],a$resid[a$color=='black'],alternative="less")$p.value, digits = 2, format = 'g'), sep='')), 
text.col = c("black","blue","red","black"))

#Fig3C
halflife = read.delim("Schofield_K562_half_lives.txt") #from Schofield et al
halflife = halflife[,c(1,6)]
colnames(halflife)[2]='half_life'
halflife$half_life=log2(halflife$half_life)
halflife = merge(halflife, c, by.x=1, by.y=2)
a1=merge(halflife, a, by.x=3, by.y=1)
"K562 half lives measured for this many genes:"
nrow(a1)
a1$quintile <- cut(a1$half_life, breaks=quantile(a1$half_life, probs=seq(0,1, by=0.2), na.rm=TRUE), include.lowest=TRUE)
boxplot(a1$resid~a1$quintile,outline=F, cex=1.5, cex.axis=2, cex.lab=2, cex.main=2, las=2, notch=T, col="red")
cor.test(a1$half_life, a1$resid)
cor.test(a1$half_life, a1$resid, method='spearman')

a$color="black"
#many groups
a$color[grep("hemoglobin subunit", a$Description)]='red'
#Xie et al 2017
a$color[a$geneName=="PIM1"]='red'
a$color[a$geneName=="SMYD3"]='red'
a$color[a$geneName=="FADS1"]='red'
a$color[a$geneName=="PRKAR2B"]='red'
#Fulco et al 2016
a$color[a$geneName=="GATA1"]='red'
a$color[a$geneName=="MYC"]='red'

#FigS3A
plot.ecdf(a$resid[a$color=='black'], xlim=c(-4,4), ylim=c(0,1), verticals= TRUE, do.points = FALSE, col="black")
plot.ecdf(a$resid[a$color=='red'], verticals= TRUE, do.points = FALSE, col="red", add=T)
legend("topleft", bg="white", bty="n", legend = c(paste("non-enhancer-driven genes, n = ", length(a$resid[a$color=='black']), sep=''),
paste("enhancer-driven genes, n = ", length(a$resid[a$color=='red']), sep=''),
paste("P value: ", formatC(ks.test(a$resid[a$color=='red'],a$resid[a$color=='black'],alternative="less")$p.value, digits = 2, format = 'g'), sep='')), text.col = c("black","red", "black"))

a$K562PredAdj=predict(lm(E123~K562Pred, a))

#Fig3A
smoothScatter(a$K562PredAdj, a$E123, cex.axis=2, cex.lab=2, bty="n", xlab="Predicted K562 expression level", ylab="K562 expression level (log10)", xlim=c(-1.5, 2), ylim=c(-1, 4), las=1, cex=.5)
abline(0,1, col="red")
text(a[a$color=="red","K562PredAdj"], a[a$color=="red","E123"], labels = a[a$color=="red","geneName"], offset = 0.5, col="red")
text(1.5, 4, labels = paste("r^2 =", round(cor(a$K562Pred, a$E123)^2,2)), offset = 0.5, col="black")

writefile(cbind(a$gene_id, a$geneName, a$Description, round(a$Pred,3), round(a$K562Pred,3), round(a$GM12878Pred,3), round(a$SuRE,3), round(a[,2:58],3)), "TableS1_human.txt")