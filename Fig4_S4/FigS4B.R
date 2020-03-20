library(data.table)
library(Biostrings)
library(rhdf5)
library(reshape2)
source("coefplot.r")

h5dir = args[1]
testIDs = h5read(paste(h5dir, "test.h5", sep='/'),"geneName")
trainIDs = h5read(paste(h5dir, "train.h5", sep='/'),"geneName")
valIDs = h5read(paste(h5dir, "valid.h5", sep='/'),"geneName")
file1 = "Roadmap_FantomAnnotations.InputData.pM10Kb.txt.gz"
file2 = "promoters_pM1.5Kb.FIMO_scanned.condensed.txt.gz"
kmerlen = 1

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

train = inp.tbl[rownames(inp.tbl) %in% trainIDs, ]
valid = inp.tbl[rownames(inp.tbl) %in% valIDs, ]
test = inp.tbl[rownames(inp.tbl) %in% testIDs, ]

mod1 = lm(EXPRESSION ~ ., data=train)
mod2 = lm(EXPRESSION ~ UTR5LEN+CDSLEN+INTRONLEN+UTR3LEN+UTR5GC+CDSGC+UTR3GC+ORFEXONDENSITY, data=train)
mod3 = lm(EXPRESSION ~ .-(UTR5LEN+CDSLEN+INTRONLEN+UTR3LEN+UTR5GC+CDSGC+UTR3GC+ORFEXONDENSITY), data=train)

cor(test$EXPRESSION, predict(mod1, newdata = test))^2
cor(test$EXPRESSION, predict(mod2, newdata = test))^2
cor(test$EXPRESSION, predict(mod3, newdata = test))^2

summary(mod2)
pdf("FigS4B_human.pdf") #change to mouse if dir is mouse
coefplot(mod2, parm = -1)
