a=read.delim(gzfile("57epigenomes.RPKM.pc.gz"))
a$median=apply(a[,3:ncol(a)],1,median)
write.table(a[,c("gene_id","median")],quote=F,row.names=F,col.names=F,sep="\t", file="57epigenomes.median_expr.txt")