library(zoo)
library(GenomicRanges)
library(mixtools)

file = args[1]
fold = args[2]
species = args[3]

if(species == 'human') expr=read.delim("all_crossvalidated_predictions.txt") else expr=read.delim("all_crossvalidated_predictions_mouse.txt")

kmerlen = 1

a=fastread(paste("zcat ", file, sep=''))
ids=a[,1]
preds = expr[match(ids, expr$Gene),"Pred"]
mixmdl = normalmixEM(preds)
thresh2 = mixmdl$mu[2]
post.df <- as.data.frame(cbind(x = mixmdl$x, mixmdl$posterior))
threshold = post.df[which(post.df$comp.1 == max(post.df$comp.1[post.df$comp.1 <= 0.5])),"x"]

# pdf("Fig6B.pdf",width=10,height=4) #general plot to make 6B histogram, but on full data rather than 1 fold
# plot(mixmdl,which=2)
# abline(v=threshold)
# dev.off()

a=a[,3:(ncol(a)-6)]

# low vs high
zeromean=apply(a[preds < threshold,], 2, mean )
zerosd=apply(a[preds < threshold,], 2, function(x) qnorm(0.995)*sd(x) ) #99th% confidence interval of z-distribution

gr <- GRanges()
gr = suppressWarnings(do.call("c", sapply(1:nrow(a), function(seq) {
    idx = which(as.vector(unlist(lapply(1:ncol(a), function(x) { a[seq,x] > zeromean[x] + zerosd[x] })))) #| a[seq,x] < zeromean[x] - zerosd[x]
    GRanges(seqnames = rep(ids[seq],length(idx)),IRanges(start = idx, end = idx+kmerlen-1))
})))

writefile(as.data.frame(reduce(gr)), paste("motif_analysis/bestpos1mer", fold, species, ".txt", sep=''), col.names=F)
