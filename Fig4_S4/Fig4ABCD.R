library(ROCR)

#### MOUSE ######

a=read.delim("all_crossvalidated_predictions_mESC.txt")
b=read.delim("all_crossvalidated_predictions_mouse.txt")
colnames(a)[2:3]=c("mESCPred","mESCActual")
a=merge(a,b,by=1)
a[a=='']="NA"
mouse = a

### HUMAN #####

a=read.delim(gzfile("57epigenomes.RPKM.pc.gz"))
a$E000=NULL
a[,2:ncol(a)]=log10(a[,2:ncol(a)]+0.1)

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

pdf("Fig4ABCD.pdf", width=8, height=8)
par(mar=c(7,7,5,5), mgp = c(5, 1.5, 0))

a$SuREAdj=predict(lm(E123~SuRE, a))
a$K562PredAdj=predict(lm(E123~K562Pred, a))
a$PromoterActivity=predict(lm(E123~SuRE+K562Pred, a))

#Fig4A
(cors = data.frame( K562=c(cor(a$K562Actual,a$Pred)^2, cor(a$K562Actual,a$K562Pred)^2),
                    GM12878=c(cor(a$GM12878Actual,a$Pred)^2, cor(a$GM12878Actual,a$K562Pred)^2),
                    mESC=c(cor(mouse$mESCActual,mouse$Pred)^2, cor(mouse$mESCActual,mouse$mESCPred)^2) ))
barplot(as.matrix(cors), las=2, beside=TRUE, col=c("red","blue"), border=F, ylim=c(0, 0.6), ylab="r^2 to gene expression level" )

#Fig4D
smoothScatter(a$SuREAdj, a$E123, cex.axis=2, cex.lab=2, bty="n", xlab="SuRE activity", ylab="K562 expression level (log10)", xlim=c(-1.5, 2), ylim=c(-1, 4), las=1, cex=.5)
abline(0,1, col="red")
text(1.5, 4, labels = paste("r^2 =", round(cor(a$E123, a$SuRE)^2,2)), offset = 0.5, col="black")

smoothScatter(a$PromoterActivity, a$E123, cex.axis=2, cex.lab=2, bty="n", xlab="Predicted expression level, joint model", ylab="K562 expression level (log10)", xlim=c(-1.5, 2), ylim=c(-1, 4), las=1, cex=.5)
abline(0,1, col="red")
text(1.5, 4, labels = paste("r^2 =", round(cor(a$PromoterActivity, a$E123)^2,2)), offset = 0.5, col="black")

#Fig4BC
a$deltaActual = a$K562Actual-a$GM12878Actual
a$deltaPred = a$K562Pred-a$GM12878Pred
nrow(a)
smoothScatter(a$GM12878Actual, a$K562Actual, cex.axis=2, cex.lab=2, bty="n", xlab="GM12878 expression level (log10)", ylab="K562 expression level (log10)", xlim=c(-1, 4), ylim=c(-1, 4), las=1, cex=.5)
abline(1,1, col="red")
abline(-1,1, col="red")
text(0, 4, labels = paste("Upregulated in K562:", nrow(a[a$deltaActual > 1,])), offset = 0.5, col="black")
text(3, -1, labels = paste("Upregulated in GM12878:", nrow(a[a$deltaActual < -1,])), offset = 0.5, col="black")

b=a[abs(a$deltaActual) > 1,]
b$deltaActual = ifelse(b$deltaActual > 0, 1, 0)
plot(performance( prediction( b$deltaPred, b$deltaActual), "tpr", "fpr"), col="blue", las=1, cex.axis=2, cex.lab=2)
text(0.2, 1, labels = paste("AUC = ", round(performance( prediction(b$deltaPred, b$deltaActual), "auc")@y.values[[1]],2), ' (n = ', nrow(b), ')', sep=''), offset = 1.5, col="black")
abline(0,1,col="grey")
dev.off()