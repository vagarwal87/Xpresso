library(data.table)
library(latticeExtra)

getbest <- function(dir){
	files = list.files(path=dir, pattern='.txt', full.names=T)
	sbtable <- as.data.frame(do.call("rbind", lapply(files, FUN=function(file){
		cmd = paste("tail -2", file, "| perl -ne \'@a=split /= /; print $a[1];\'")
		tmp = t(read.table(textConnection(system(cmd, intern=TRUE))))
		tmp$file = strsplit(file, "/")[[1]][2]
		names(tmp) = c("r2","MSE","samples")
		tmp
	})))

	sbtable$samples = as.character(sbtable$samples)
	sbtable = data.table(sbtable)
	sbtable = sbtable[ , .SD[which.min(MSE)], by = samples]
	sbtable
}

a=list()
a[[1]]=getbest("ortholog_results/train_human_test_human/")
a[[2]]=getbest("ortholog_results/train_human_test_mouse/")
a[[3]]=getbest("ortholog_results/train_mouse_test_human/")
a[[4]]=getbest("ortholog_results/train_mouse_test_mouse/")

a=as.data.frame(do.call("rbind", a))
a$samples=as.factor(a$samples)

a
pdf("Fig2D.pdf", height=4, width=5)
obj1 = xyplot(MSE ~ samples, a, type = "p", pch=19, lwd=2, scales=list(x=list(rot=45)))
obj2 = xyplot(r2 ~ samples, a, type = "p", pch=19, lwd=2, scales=list(x=list(rot=45)))
doubleYScale(obj1, obj2, add.ylab2 = TRUE)
dev.off()
