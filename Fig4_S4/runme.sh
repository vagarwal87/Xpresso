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

#generate baseline dinucleotide model, and prepare to extract features for PWM-based model
zcat Roadmap_FantomAnnotations.InputData.pM10Kb.txt.gz | cut -f 1,11 | tail -n+2 | perl -ne '@a=split; print ">$a[0]\n".substr($a[1],8500,3000)."\n";' >promoters_pM1.5Kb.fa
zcat Mouse_FantomAnnotations.InputData.pM10Kb.txt.gz | cut -f 1,11 | tail -n+2 | perl -ne '@a=split; print ">$a[0]\n".substr($a[1],8500,3000)."\n";' >promoters_pM1.5Kb.mouse.fa
fasta-get-markov promoters_pM1.5Kb.fa >promoters_pM1.5Kb.firstOrderMarkov_background
fasta-get-markov promoters_pM1.5Kb.mouse.fa >promoters_pM1.5Kb.mouse.firstOrderMarkov_background

fimo --bgfile promoters_pM1.5Kb.firstOrderMarkov_background --verbosity 1 --text --skip-matched-sequence JASPAR_CORE_2016_vertebrates.meme promoters_pM1.5Kb.fa | gzip -c >promoters_pM1.5Kb.FIMO_scanned.txt.gz
fimo --bgfile promoters_pM1.5Kb.mouse.firstOrderMarkov_background --verbosity 1 --text --skip-matched-sequence JASPAR_CORE_2016_vertebrates.meme promoters_pM1.5Kb.mouse.fa | gzip -c >promoters_pM1.5Kb.mouse.FIMO_scanned.txt.gz

zcat promoters_pM1.5Kb.FIMO_scanned.txt.gz | cut -f 1,2 | uniq | gzip -c >promoters_pM1.5Kb.FIMO_scanned.condensed.txt.gz
zcat promoters_pM1.5Kb.mouse.FIMO_scanned.txt.gz | cut -f 1,2 | uniq | gzip -c >promoters_pM1.5Kb.mouse.FIMO_scanned.condensed.txt.gz

#generate 1mer-6mer baseline models for human and mouse, respectively
for x in {1..6}; do { Rscript baseline_models.R pM10Kb_1KTest $x; } done &
for x in {1..6}; do { Rscript baseline_models.R pM10Kb_1KTest_Mouse $x; } done &

# empirical results from Xpresso and baselines are stored in model_comparison.txt
Rscript Fig4E_S4.R

#generates for human, can change directory and code to generate for mouse
Rscript FigS4B.R

#generates for human, can change directory and code to generate for mouse
Rscript FigS4C.R
