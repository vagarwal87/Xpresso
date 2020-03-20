library(data.table)
library(latticeExtra)

setwd("subsampled_10fold/")

files = list.files(path='.', pattern='trial', recursive=T)
# files = files[grepl("^2000|4000|6000|8000|10000|14000|16000",files)]
table <- as.data.frame(do.call("rbind", lapply(files, FUN=function(file){
	cmd = paste("tail -2", file, "| perl -ne \'@a=split /= /; print $a[1];\'")
	tmp = t(read.table(textConnection(system(cmd, intern=TRUE))))
	tmp$sample = as.numeric(dirname(file))
	tmp$rep = as.numeric(strsplit(basename(file), "_")[[1]][1])
	tmp$trial = as.numeric(strsplit(strsplit(file, "_trial")[[1]][2], '\\.')[[1]][1])
	names(tmp) = c("r2","MSE","sample","rep","trial")
	tmp
})))

table = as.data.frame(apply(table,2,function(x) as.numeric(as.character(x))))
table = do.call("rbind", lapply(unique(table$sample), function(sub) { do.call("rbind", lapply(unique(table$rep), function(x) { tmp=table[table$rep==x & table$sample==sub,]; tmp[which( tmp$MSE == min(tmp$MSE) ),] } )) }) )

table = as.data.frame(aggregate(.~sample,table,function(x) c(mean=mean(x), sd=sd(x))))
table
table[,2][,2]=table[,2][,2]/sqrt(10) #std err
table[,3][,2]=table[,3][,2]/sqrt(10) #std err
table

pdf("subsample.pdf", height=3, width=5)
obj1 = xyplot(MSE[,1] ~ sample, table,
       panel = function(x, y, ...){
         panel.arrows(x, y, x, table[,3][,1]+1.96*table[,3][,2], length = 0, angle = 90)
         panel.arrows(x, y, x, table[,3][,1]-1.96*table[,3][,2], length = 0, angle = 90)
         panel.xyplot(x, y, ...)
	 }, type = "o" , ylim=c(0.4,0.5), lwd=2, scales = list(x = list(at = seq(2000,16000,2000) ))) #limits = c(-500,500)
obj2 = xyplot(r2[,1] ~ sample, table,
       panel = function(x, y, ...){
         panel.arrows(x, y, x, table[,2][,1]+1.96*table[,2][,2], length = 0, angle = 90)
         panel.arrows(x, y, x, table[,2][,1]-1.96*table[,2][,2], length = 0, angle = 90)
         panel.xyplot(x, y, ...)
	 }, type = "o", lwd=2, ylim=c(0.5,0.6))
doubleYScale(obj1, obj2, add.ylab2 = TRUE)
dev.off()
