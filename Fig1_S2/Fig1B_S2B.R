library(latticeExtra)

getresults = function(thisfile){
    sites = read.table(text = system(paste("python print_losses.py", thisfile), intern=T), sep='\t')
    colnames(sites)=c("leftpos","rightpos","loss","params")
    print(sites[which(sites$loss == min(sites$loss)),])
    say(nrow(sites), "trials")
    sites
}

c = getresults(args[1])
e = getresults(args[2])

pdf("Fig1B.pdf",width=5,height=4)
plot(1:nrow(c), sapply(1:nrow(c), function(x) min(c[1:x, "loss"])), lwd=2, bty='n', col='red', type="l", 
xlim = c(0, 1000), ylim = c(0.4, 0.7), xlab="Number of iterations", ylab="Validation mean squared error, best model found")
abline(h=0.479, lwd=2, lty=2, col='black')
lines(1:nrow(e), sapply(1:nrow(e), function(x) min(e[1:x, "loss"])), lwd=2, col='purple')
legend("topright", bg="white", bty="n", legend = c("Tree of Parzen estimators", "Simulated annealing","Best manually discovered, -1.5Kb to 1.5Kb"),
text.col = c("red", "purple","black"), cex=0.8)
dev.off()

plotboundaries = function(a){
    b=aggregate(a$loss, by=list(leftpos=a$leftpos, rightpos=a$rightpos), min)
    totsize = 10000
    pdf("FigS2B.pdf",width=5,height=6)
    layout(matrix(c(1,1,1,2), 1, 4))
    par(mar = c(1, 1, 5, 1))
    b$mycol = as.character("red")
    N=min(nrow(b),100)
    b=b[order(b$x, decreasing=T),]
    b=b[(nrow(b)-N+1):nrow(b),]
    b=rbind(c(8500,11500,0.479,"blue"), b)
    b$leftpos=as.integer(b$leftpos) - 10000
    b$rightpos=as.integer(b$rightpos) - 10000
    plot(1:(N+1),xlim=c(-totsize,totsize), type="n", cex.lab = 2, bty="n", yaxt='n', xaxt='n')
    axis(3, at=seq(-totsize,totsize,totsize/5))
    mtext("Position relative to TSS", side=3, line=3)
    for(x in 1:(N+1)) lines(c(b$leftpos[x],b$rightpos[x]), c(x,x), col=b$mycol[x], type="l", lty=1, lwd=2)
    abline(v=0, lwd=2, col='black')

    plot(1:(N+1),xlim=c(0.4,0.48), type="n", cex.lab = 2, bty="n", yaxt='n', xaxt='n')
    axis(3, at=seq(0.4,0.48,0.02))
    mtext("Validation MSE", side=3, line=3)
    for(x in 1:(N+1)) lines(c(0,b$x[x]), c(x,x), col="grey", type="l", lty=1, lwd=2)
    abline(v=0.479, lwd=2, lty=2, col='blue')
    dev.off()
}

plotboundaries(c)