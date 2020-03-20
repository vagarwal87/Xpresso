library(latticeExtra)

crp.rg <- colorRampPalette(c("red","orange","green","cyan","blue","purple","magenta"))
cols <- sample(crp.rg(10))

plotresults = function(dir){
    files = paste(list.files(path=dir, pattern='.txt', full.names=T))
    pdf("Fig1C.pdf",width=5,height=4)
    plot(0, lwd=1, bty='n', type="l", xlim = c(0, 25), ylim = c(0.4, 1), xlab="Epoch", ylab="Validation MSE", las = 1)
    abline(h=0.479, lwd=1, lty=2, col='black')
    lapply(1:length(files), FUN=function(i){
        file = files[i]
        cmd = paste("grep val_loss", file, " | perl -ne 'chomp; ($mse) = ($_ =~ /val_loss: (\\d+.\\d+)/); print \"$mse \";'")
        sites = unlist(strsplit(system(cmd, intern=T), "\\s+"))
        sites = as.numeric(sites[2:length(sites)])
        lines(1:length(sites), sites, lwd=2, col=cols[i])
    })
    dev.off()
}

c = plotresults(args[1])