library(data.table)
library(Biostrings)
library(rhdf5)
library(reshape2)
library(beeswarm)

h5dir = args[1]

file1 = "Roadmap_FantomAnnotations.InputData.pM10Kb.txt.gz"
file2 = "promoters_pM1.5Kb.FIMO_scanned.condensed.txt.gz"
kmerlen = 4

if (grepl("Mouse",h5dir)) {
    file1 = "Mouse_FantomAnnotations.InputData.pM10Kb.txt.gz"
    file2 = "promoters_pM1.5Kb.mouse.FIMO_scanned.condensed.txt.gz"
}

inp.tbl <- fread(paste("zcat", file1),header=T,data.table=F,sep="\t")
rownames(inp.tbl) = inp.tbl[,1]
inp.tbl[,1] = NULL

inp.tbl$PROMOTER = substring(inp.tbl$PROMOTER,8500,11500)

inp.tbl=cbind(inp.tbl, do.call(rbind, lapply(inp.tbl$PROMOTER, function(x){
    y = oligonucleotideFrequency(DNAStringSet(x), kmerlen)
    y/sum(y)
})))
inp.tbl$T=NULL #remove TT dinucleotide to ensure matrix is full rank

inp.tbl$PROMOTER = NULL
inp.tbl[,c(1:5, 9)] = log10(inp.tbl[,c(1:5, 9)]+0.1)
inp.tbl=as.data.frame(scale(inp.tbl))

# save(inp.tbl, file="5merInputTable.RData")
# # load("TriInputTable.RData")

motif_hits <- fread(paste("zcat", file2),header=T,data.table=F,sep="\t")
motif_hits <- dcast(motif_hits, motif_hits[,2] ~ motif_hits[,1], function(x) 1, fill = 0)
inp.tbl=merge(inp.tbl, motif_hits, by.x=0, by.y=1, all.x=T)

sum(is.na(inp.tbl))
inp.tbl[is.na(inp.tbl)] = 0
rownames(inp.tbl) = inp.tbl[,1]
inp.tbl[,1] = NULL

z = do.call("rbind", lapply(1:10, function(i){
    testIDs = h5read(paste(h5dir, '/', i, "test.h5", sep=''),"geneName")
    trainIDs = h5read(paste(h5dir, '/', i, "train.h5", sep=''),"geneName")
    valIDs = h5read(paste(h5dir, '/', i, "valid.h5", sep=''),"geneName")
    train = inp.tbl[rownames(inp.tbl) %in% trainIDs | rownames(inp.tbl) %in% valIDs, ]
    # valid = inp.tbl[rownames(inp.tbl) %in% valIDs, ]
    test = inp.tbl[rownames(inp.tbl) %in% testIDs, ]

    mod1 = lm(EXPRESSION ~ ., data=train)
    c(i, cor(test$EXPRESSION, predict(mod1, newdata = test))^2)
}))
colnames(z)=c("fold","baseliner2")

cv_folder = paste0(h5dir, 2)
files = list.files(path=cv_folder, pattern='.txt.gz', full.names=T)
table <- do.call("rbind", lapply(files, FUN=function(file){
        cmd = paste("zcat ", file, "| tail -2 | perl -ne \'@a=split /= /; print $a[1];\'")
        tmp = t(read.table(textConnection(system(cmd, intern=TRUE))))
        tmp$fold = as.numeric(strsplit(basename(file), "_")[[1]][1])
        tmp$trial = as.numeric(strsplit(strsplit(file, "_trial")[[1]][2], '\\.')[[1]][1])
        names(tmp) = c("r2","MSE","fold","trial")
        tmp
}))

table=as.data.frame(apply(table,2,function(x) as.numeric(as.character(x))))
table = do.call("rbind",lapply(unique(table$fold), function(x) { tmp=table[table$fold==x,]; tmp[which( tmp$MSE==min(tmp$MSE) ),c("r2","fold")] } ))
table = aggregate(table$r2, by=list(fold=table$fold), mean)
colnames(table)[2]="Xpressor2"
table

cvr2 = merge(table,z,by=1)
t.test(cvr2[,2],cvr2[,3],paired=T)
cvr2
pdf("FigS4C_human.pdf") #change to mouse if dir is mouse
beeswarm(cvr2[,2:3],ylim=c(0,1), las=2, bty='n', pch=19) #0.4,0.8
dev.off()
