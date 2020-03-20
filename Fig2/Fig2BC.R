library(LSD)
library(data.table)

#############################
##Species-specific analysis##
#############################

pdf("Fig2B.pdf")
a=read.delim("../datasets/pM10Kb_1KTest/predictions.txt")
actual=read.delim("57epigenomes.median_expr.txt",F)
colnames(actual)=c("Gene","UnscaledExpr")
actual$UnscaledExpr=log10(actual$UnscaledExpr+0.1)
a=merge(a,actual,by=1)
model=lm(UnscaledExpr~Actual,a)
a$Pred=predict(model,newdata=data.frame(Actual=a$Pred))
a$Actual=predict(model)

"Human r^2:"
cor(a$Pred,a$Actual)^2
heatscatter(a$Pred, a$Actual, bty='n', xlim=c(-1,3), ylim=c(-1,3), cex.axis=2, cex.lab=2, las=1)
a=read.delim("../datasets/pM10Kb_1KTest_Mouse/predictions.txt")
actual=read.delim("mouse.median_expr.txt",F)
colnames(actual)=c("Gene","UnscaledExpr")
actual$UnscaledExpr=log10(actual$UnscaledExpr+0.1)
a=merge(a,actual,by=1)
model=lm(UnscaledExpr~Actual,a)
a$Pred=predict(model,newdata=data.frame(Actual=a$Pred))
a$Actual=predict(model)
dev.off()

"Mouse r^2:"
cor(a$Pred,a$Actual)^2
pdf("Fig2C.pdf")
heatscatter(a$Pred, a$Actual, bty='n', xlim=c(-1,3), ylim=c(-1,3), cex.axis=2, cex.lab=2, las=1)
dev.off()
