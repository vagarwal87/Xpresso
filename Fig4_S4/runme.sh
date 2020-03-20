
perl -ne 'print "chr$_";' hg38_promoters_cage_corrected.bed >hg38_promoters_cage_corrected_withChr.bed
liftOver -bedPlus=6 hg38_promoters_cage_corrected_withChr.bed hg38ToHg19.over.chain hg19_promoters_cage_corrected_withChr.bed unmapped
./supplement_ids.pl >hg19_promoters_cage_corrected_withChr_andOthers.bed
grep -P '\-$' hg19_promoters_cage_corrected_withChr_andOthers.bed >hg19_promoters_cage_corrected_withChr_andOthers_minus.bed
grep -P '\+$' hg19_promoters_cage_corrected_withChr_andOthers.bed >hg19_promoters_cage_corrected_withChr_andOthers_plus.bed

wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE78nnn/GSE78709/suppl/GSE78709_sure23.plasmid.norm.combined.45.55.plus.160504.bw
wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE78nnn/GSE78709/suppl/GSE78709_sure23.plasmid.norm.combined.45.55.minus.160504.bw

bigWigAverageOverBed -sampleAroundCenter=1000 GSE78709_sure23.plasmid.norm.combined.45.55.plus.160504.bw hg19_promoters_cage_corrected_withChr_andOthers_plus.bed GSE78709_sure23.plasmid.norm.combined.45.55.plus.promoters.bigWigSignal
bigWigAverageOverBed -sampleAroundCenter=1000 GSE78709_sure23.plasmid.norm.combined.45.55.minus.160504.bw hg19_promoters_cage_corrected_withChr_andOthers_minus.bed GSE78709_sure23.plasmid.norm.combined.45.55.minus.promoters.bigWigSignal

Rscript Fig4ABCD.R

Rscript Fig4E_S4.R

#generates for human, can change directory and code to generate for mouse
Rscript FigS4B.R

#generates for human, can change directory and code to generate for mouse
Rscript FigS4C.R
