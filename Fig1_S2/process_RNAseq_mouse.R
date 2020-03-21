a=read.delim(gzfile("mouse_FPKMs.tsv.gz"), F)
colnames(a)[1]="gene_id"
a[,1]=substring(a[,1],1,18)
a$median=apply(a[,2:ncol(a)],1,median)
write.table(a[,c("gene_id","median")],quote=F,row.names=F,col.names=F,sep="\t", file="mouse.median_expr.txt")