options(warn=-1)

a=read.delim("model_comparison.txt", sep=' ')
a

pdf("FigS4A.pdf", height=8, width=10)
par(oma=c(1,20,1,1))
barplot(rbind(a$test_r_squared,a$test_r_squared_withHL),beside=T,horiz=T,
names.arg=a$model,las=1,col=c("red","blue"), border=F, xlim=c(0,0.8))
legend("bottomright", bg="white", bty="n", legend = c("with half life", "without half life"), text.col = c("blue","red"))
dev.off()

a=read.delim("model_comparison_Fig3.txt", sep=' ')
a
pdf("Fig4E.pdf", height=8, width=10)
par(oma=c(1,20,1,1))
barplot(a$r_squared,horiz=T, names.arg=a$model,las=1,col=c(rep("red",11),rep("blue",22)), border=F, xlim=c(0,0.8))
legend("bottomright", bg="white", bty="n", legend = c("mouse", "human"), text.col = c("red","blue"))
dev.off()
