library(data.table)
library(latticeExtra)

cv_folder = args[1]
predfolder = args[2]
outfile = args[3]

files = list.files(path=cv_folder, pattern='.txt', full.names=T)

table <- do.call("rbind", lapply(files, FUN=function(file){
	cmd = paste("tail -2", file, "| perl -ne \'@a=split /= /; print $a[1];\'")
	tmp = t(read.table(textConnection(system(cmd, intern=TRUE))))
	tmp$fold = as.numeric(strsplit(basename(file), "_")[[1]][1])
	tmp$trial = as.numeric(strsplit(strsplit(file, "_trial")[[1]][2], '\\.')[[1]][1])
	names(tmp) = c("r2","MSE","fold","trial")
	tmp
}))

table=as.data.frame(apply(table,2,function(x) as.numeric(as.character(x))))
head(table)

do.call("rbind", lapply(unique(table$fold), function(x) { tmp=table[table$fold==x,]; tmp[which( tmp$MSE == min(tmp$MSE) ),] } ) )
table = do.call("rbind",lapply(unique(table$fold), function(x) { tmp=table[table$fold==x,]; tmp[which( tmp$MSE==min(tmp$MSE) ),c("fold","trial")] } ))

if (nrow(table) == 10){
    files = apply(table, 1, function(x) { paste(predfolder,x[2],x[1],"predictions.txt",sep='') } )
    say(files)
    table = do.call("rbind", lapply(files, function(x) { read.delim(x) } ) )
    write.table(table,file=outfile, quote=F, row.names=F, sep='\t')
}
#otherwise cant do, select which trial to use from table due to tie
