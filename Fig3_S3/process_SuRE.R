library(data.table)
a=fread("zcat GSE78709_SuRE-counts-K562_B45_B55_LP170105.txt.gz",data.table=F,sep='\t')
a$rep1cDNA= (a$"SuRE K562 cDNA, Biorep1, PCRrep1"+a$"SuRE K562 cDNA, Biorep1, PCRrep2"+a$"SuRE K562 cDNA, Biorep1, PCRrep3"+a$"SuRE K562 cDNA, Biorep1, PCRrep4")
a$rep2cDNA= (a$"SuRE K562 cDNA, Biorep2, PCRrep1"+a$"SuRE K562 cDNA, Biorep2, PCRrep2")
a$seqnum=paste("seq",1:nrow(a),sep='')
a$empty='.'
a=a[,c("chr","start","end","seqnum","empty","strand","rep1cDNA","rep2cDNA","plasmid library, PCR rep1","plasmid library, PCR rep2")]
writefile(a,"GSE78709_SuRE-counts-K562_B45_B55_LP170105.bed",col.names=F)