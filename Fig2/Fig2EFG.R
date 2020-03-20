library(LSD)
library(data.table)
library(ROCR)

####################
##One2One analysis##
####################
a=read.delim("human2mouse_one2one_orthologs.txt",F)
b=fread("zcat Roadmap_FantomAnnotations.InputData.pM10Kb.txt.gz | cut -f 1-2",header=T,data.table=F,sep="\t")
c=fread("zcat Mouse_FantomAnnotations.InputData.pM10Kb.txt.gz | cut -f 1-2",header=T,data.table=F,sep="\t")
b$EXPRESSION = log10(b$EXPRESSION+0.1)
c$EXPRESSION = log10(c$EXPRESSION+0.1)
d=merge(a,c,by.x=2,by.y=1)
d=merge(d,b,by.x=2,by.y=1)

colnames(d)=c("hid","mid","type","mexpr","hexpr")
head(d)

say(nrow(d), "genes")
say("Pearson corr =", cor(d$hexpr,d$mexpr))

pdf("Fig2E.pdf")
par(mar=c(7,7,5,5), mgp = c(5, 1, 0))
plot.ecdf(b$EXPRESSION, xlim=c(-1,3), ylim=c(0,1), verticals= TRUE, do.points = FALSE, col="purple",
    ylab="Cumulative fraction", xlab="log10(Median expression level + 0.1)", col.01line = "white", lwd=2, cex.axis=2, cex.lab=2, bty="n", las=1)
plot.ecdf(c$EXPRESSION, verticals= TRUE, do.points = FALSE, add = TRUE, col = "red", col.01line = "white", lwd=2, bty="n")
plot.ecdf(d$hexpr, verticals= TRUE, do.points = FALSE, add = TRUE, col = "blue", col.01line = "white", lwd=2, bty="n")
plot.ecdf(d$mexpr, verticals= TRUE, do.points = FALSE, add = TRUE, col = "cyan", col.01line = "white", lwd=2, bty="n")

legend("bottomright", bg="white", bty="n", legend =
  c( paste("human (", length(b$EXPRESSION), ")") , paste("mouse (", length(c$EXPRESSION), ")"),
  paste("human, one-to-one orthologs (", length(d$hexpr), ")"), paste("mouse, one-to-one orthologs (", length(d$mexpr), ")")),
  text.col = c("purple","red","blue","cyan"))
dev.off()

pdf("Fig2F.pdf")
heatscatter(d$hexpr, d$mexpr, xlab="Human", ylab="Mouse", bty='n', cex=0.3, xlim=c(-1,3), ylim=c(-1,3), cex.axis=2, cex.lab=2, las=1)
dev.off()

writefile(d,"1to1_orthologs_expression.txt", col.names=F)

pdf("Fig2G.pdf")
b=read.delim("all_crossvalidated_predictions.txt")
c=read.delim("all_crossvalidated_predictions_mouse.txt")
d=merge(d,b,by=1)
colnames(d)[4:5]=c("mouse_expr","human_expr")
colnames(d)[6:7]=c("human_pred","human_Actual")
e=merge(d,c,by.x=2,by.y=1)
colnames(e)[8:9]=c("mouse_pred","mouse_Actual")
attach(e)
head(e)

e$diff = mouse_expr-human_expr
f=e[abs(e$diff) > 1,]
f$mouseOrHuman = ifelse(f$diff > 0, 1, 0)
"human-specific"
sum(f$diff < 0)
"mouse-specific"
sum(f$diff > 0)
head(f)
plot(performance( prediction( f$mouse_pred-f$human_pred, f$mouseOrHuman), "tpr", "fpr"), col="blue", las=1, cex.axis=2, cex.lab=2, bty='n')
text(0.2, 1, labels = paste("AUC = ", round(performance( prediction(f$mouse_pred-f$human_pred, f$mouseOrHuman), "auc")@y.values[[1]],2), ' (n = ', nrow(f), ')', sep=''), offset = 1.5, col="black")
abline(0,1,col="grey")
dev.off()
