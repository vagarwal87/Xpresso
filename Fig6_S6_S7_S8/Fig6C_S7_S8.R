file = args[1]
outfile = args[2]

b=list()
mononuc=list()

################ cpg
a=fastread(paste("zcat ", file,sep=''))
a$PROMOTER = substring(a$PROMOTER,6000,14000)
a$EXPRESSION=log10(a$EXPRESSION+0.1)
a$EXPRESSION=(a$EXPRESSION-min(a$EXPRESSION))/(max(a$EXPRESSION)-min(a$EXPRESSION))
a$bin <- cut(a$EXPRESSION, breaks=c(0,0.001,0.33,0.66,1), include.lowest = TRUE)

# b=lapply(levels(a$bin), function(x){ table(do.call(c, lapply(a[a$bin==x,"PROMOTER"], function(y) unlist(gregexpr("CG|cG|Gc|cg",y)) )))/nrow(a[a$bin==x,]) } )
# g=lapply(levels(a$bin), function(x){ table(do.call(c, lapply(a[a$bin==x,"PROMOTER"], function(y) unlist(gregexpr("G|g",y)) )))/nrow(a[a$bin==x,]) } )
# c=lapply(levels(a$bin), function(x){ table(do.call(c, lapply(a[a$bin==x,"PROMOTER"], function(y) unlist(gregexpr("C|c",y)) )))/nrow(a[a$bin==x,]) } )

for (nuc1 in c('a','c','g','t')){
	say(nuc1)
	mononuc[[nuc1]]=lapply(levels(a$bin), function(x){ table(do.call(c, lapply(a[a$bin==x,"PROMOTER"], function(y) unlist(gregexpr(paste(nuc1,'|',toupper(nuc1),sep=''),y)) )))/nrow(a[a$bin==x,]) } )
}

pdf(outfile,width=10,height=8) #makes Fig6C/FigS7 for human for FigS8 for mouse
par(mfrow=c(4,4), oma = c(5,4,0,0) + 0.1, mar = c(0,0,1,1) + 0.1)
i=0
for (nuc1 in c('a','c','g','t')){
	for (nuc2 in c('a','c','g','t')){
		dinuc = paste(nuc1,nuc2,'|',toupper(nuc1),toupper(nuc2),'|',nuc1,toupper(nuc2),'|',toupper(nuc1),nuc2,sep='')
		dinuc2 = paste(toupper(nuc1),toupper(nuc2),sep='')
		say(dinuc)
		b=lapply(levels(a$bin), function(x){ table(do.call(c, lapply(a[a$bin==x,"PROMOTER"], function(y) unlist(gregexpr(dinuc,y)) )))/nrow(a[a$bin==x,]) } )

		x = 1:8000
		xin = -4000:3999
		idx = 1:8000

		plot(xin, predict(loess(b[[4]][idx]/(mononuc[[nuc1]][[4]][idx]*mononuc[[nuc2]][[4]][idx+1])~x, span=0.01), newdata=idx), col='cyan', type='l', ylim=c(0,1.6), main = dinuc2, axes = FALSE)
		axis(side = 1, labels = (i %/% 4 == 3))
		axis(side = 2, labels = (i %% 4 == 0), las=1)
		lines(xin, predict(loess(b[[3]][idx]/(mononuc[[nuc1]][[3]][idx]*mononuc[[nuc2]][[3]][idx+1])~x, span=0.01), newdata=idx), col='blue', type='l')
		lines(xin, predict(loess(b[[2]][idx]/(mononuc[[nuc1]][[2]][idx]*mononuc[[nuc2]][[2]][idx+1])~x, span=0.01), newdata=idx), col='red', type='l')
		lines(xin, predict(loess(b[[1]][idx]/(mononuc[[nuc1]][[1]][idx]*mononuc[[nuc2]][[1]][idx+1])~x, span=0.01), newdata=idx), col='black', type='l')
		i=i+1
	}
}
title(xlab = "Position relative to TSS", ylab = "Observed/Expected", outer = TRUE, line = 3)
# legend("topleft", bg="white", bty="n", legend = levels(a$bin), text.col = c("black", "red","blue","cyan"), outer = TRUE)
dev.off()
