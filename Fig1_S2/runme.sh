########### MOST OF THESE STEPS HAVE PRECOMPUTED RESULTS

#run the following in the base Xpresso folder to retrieve these precomputed results:
wget -r -np -nH --reject "index.html*" --cut-dirs 5 https://krishna.gs.washington.edu/content/members/vagar/Xpresso/data/datasets/

########### EXTRACT PROMOTERS FROM FANTOM5 CAGE PEAKS ############

# download human CAGE annotations
# the original link (now broken, stored in "datasets/ was downloaded from here)
# wget http://fantom.gsc.riken.jp/5/datahub/hg38/peaks/hg38.cage_peak_phase1and2combined.bb http://fantom.gsc.riken.jp/5/datahub/mm10/peaks/mm10.cage_peak_phase1and2combined.bb

# the revised link can be found here:
wget http://fantom.gsc.riken.jp/5/datahub/hg38/peaks/hg38.cage_peak.bb http://fantom.gsc.riken.jp/5/datahub/mm10/peaks/mm10.cage_peak.bb
bigBedToBed hg38.cage_peak_phase1and2combined.bb hg38.cage_peak_phase1and2combined.bed
bigBedToBed mm10.cage_peak_phase1and2combined.bb mm10.cage_peak_phase1and2combined.bed

# extracts best peak for each gene, removes chrY/M genes
grep -e 'p1@' hg38.cage_peak_phase1and2combined.bed | \
    perl -ne 'chomp; @a=split /\t/; $a[0]=substr($a[0],3); $mid = int($a[-2]); $start=$mid-10000; $stop=$mid+10000; print join("\t",$a[0],$start,$stop,$a[3],0,$a[5])."\n" if $start > 0;' | \
    grep -v -P "^Y|^M" >hg38_cage_promoters.bed
grep -e 'p1@' mm10.cage_peak_phase1and2combined.bed | \
    perl -ne 'chomp; @a=split /\t/; $mid = int($a[-2]); $start=$mid-10000; $stop=$mid+10000; print join("\t",$a[0],$start,$stop,$a[3],0,$a[5])."\n" if $start > 0;' | \
    grep -v -P "^chrY|^chrM" >mm10_cage_promoters.bed

#acquired 2 additional tables from BioMart and HGNC in addition to these for the mouse
wget http://www.informatics.jax.org/downloads/reports/MGI_Gene_Model_Coord.rpt http://www.informatics.jax.org/downloads/reports/MGI_EntrezGene.rpt
# converts IDs of protein-coding genes into Ensembl IDs for top promoter CAGE peak
./geneName2Ensembl.pl
# acquire Ouyang_mESC_RPKM.txt from Supplementary table of Ouyang et al.
./geneName2EnsemblMouse.pl

# extract CAGE-revised promoter sequence from hg38 genome -- REQUIRES DOWNLOAD OF HG38 AND MM10 GENOMES
bedtools getfasta -s -name -fi human_hs38_noAlt/whole_genome.fa -bed hg38_cage_promoters_ensemblID.bed -fo hg38_cage_promoters_ensemblID.fa
gzip hg38_cage_promoters_ensemblID.fa
bedtools getfasta -s -name -fi mus_musculus/mm10.fa -bed mm10_cage_promoters_ensemblID.bed -fo mm10_cage_promoters_ensemblID.fa
gzip mm10_cage_promoters_ensemblID.fa


########### EXTRACT PROMOTERS FROM ENSEMBL ############

# download human/mouse gene annotations on hg38/mm10
wget ftp://ftp.ensembl.org/pub/release-90/gtf/homo_sapiens/Homo_sapiens.GRCh38.90.gtf.gz
wget ftp://ftp.ensembl.org/pub/release-90/gtf/mus_musculus/Mus_musculus.GRCm38.90.gtf.gz
#Ensembl 90 on hg19
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_27/GRCh37_mapping/gencode.v27lift37.basic.annotation.gtf.gz

# choose 1 representative transcript for each protein-coding gene, keep chrX or chr[1..22] genes only
./choose_reference_genes.pl Homo_sapiens.GRCh38.90.gtf.gz | gzip -c >Homo_sapiens.GRCh38.90.chosenTranscript.gtf.gz
./choose_reference_genes_forhg19.pl gencode.v27lift37.basic.annotation.gtf.gz | gzip -c >Homo_sapiens.hg19.90.chosenTranscript.gtf.gz
./choose_reference_genes.pl Mus_musculus.GRCm38.90.gtf.gz | gzip -c >Mus_musculus.GRCm3.90.chosenTranscript.gtf.gz
zgrep transcript Homo_sapiens.hg19.90.chosenTranscript.gtf.gz | gzip -c >Homo_sapiens.hg19.90.chosenTranscript.geneBoundaries.gtf.gz

# generate input file for UCSC genome browser, extract 5' UTR, ORF, and 3' UTR sequences using these files using the Table Browser
perl -ne '@a = split; $a[-1] = "gene_id \"$a[-1]\"; transcript_id \"$a[-1]\""; print "chr".join("\t", @a), "\n";' \
    <(zcat Homo_sapiens.GRCh38.90.chosenTranscript.gtf.gz) | \
    gzip -c >Homo_sapiens.GRCh38.90.chr.gtf.gz
perl -ne '@a = split; $a[-1] = "gene_id \"$a[-1]\"; transcript_id \"$a[-1]\""; print "chr".join("\t", @a), "\n";' \
    <(zcat Mus_musculus.GRCm3.90.chosenTranscript.gtf.gz) | \
    gzip -c >Mus_musculus.GRCm3.90.chr.gtf.gz

# process into BED and extract +/- 10Kb region surrounding TSS
./extract_promoters.pl Homo_sapiens.GRCh38.90.chosenTranscript.gtf.gz >hg38_promoters.bed
./extract_promoters.pl Mus_musculus.GRCm3.90.chosenTranscript.gtf.gz mouse >mm10_promoters.bed

# extract Ensembl-annotated promoter sequence from hg38 genome
bedtools getfasta -s -name -fi human_hs38_noAlt/whole_genome.fa -bed hg38_promoters.bed -fo hg38_promoters.fa
gzip hg38_promoters.fa
bedtools getfasta -s -name -fi mus_musculus/mm10.fa -bed mm10_promoters.bed -fo mm10_promoters.fa
gzip mm10_promoters.fa

# histone genes to filter out, not quantified correctly due to lack of poly(A) tail
grep HIST ensembl2geneName_v90.txt | cut -f 1 >mask_histone_genes.txt
grep Hist ensembl2geneName_v90_mm10.txt | cut -f 1 >mask_histone_genes_mm10.txt

########### COLLECT & PROCESS GENE EXPRESSION DATA ############

# download pre-processed RNA-seq data from 56 cell types (+1 universal reference)
wget http://egg2.wustl.edu/roadmap/data/byDataType/rna/expression/57epigenomes.RPKM.pc.gz http://egg2.wustl.edu/roadmap/data/byDataType/rna/expression/EG.name.txt
# extract median expression values
Rscript process_RNAseq.R
cut -f 1,56 <(zcat 57epigenomes.RPKM.pc.gz) | tail -n+2 >57epigenomes.K562.txt
cut -f 1,50 <(zcat 57epigenomes.RPKM.pc.gz) | tail -n+2 >57epigenomes.GM12878.txt

cut -f 11 files_geneQuant_Rep1.txt | tail -n+2 >urls.txt
while read p; do qsub "wget $p"; done <urls.txt
gzip *.tsv
for x in *.tsv.gz; do { X=`basename $x .tsv.gz`; zgrep -P 'ENS' $x | cut -f 1,7 | gzip -c >$X.FPKM.tsv.gz; } done
paste *tsv | cut -f 1,$(echo `seq 2 2 1000` | perl -ne '@a=split / /, $_; print join(",",@a);') | sort | gzip -c >mouse_FPKMs.tsv.gz
Rscript process_RNAseq_mouse.R #generates mouse.median_expr.txt

########### GENERATE TRAINING/VALIDATION/TEST SET AND OPTIMIZE ############

#generate training/validation/test sets
perl generate_training_input.pl 57epigenomes.median_expr.txt | gzip -c >Roadmap_FantomAnnotations.InputData.pM10Kb.txt.gz
perl generate_training_input.pl 57epigenomes.K562.txt | gzip -c >Roadmap_FantomAnnotations.InputData.pM10Kb.K562expr.txt.gz
perl generate_training_input.pl 57epigenomes.GM12878.txt | gzip -c >Roadmap_FantomAnnotations.InputData.pM10Kb.GM12878expr.txt.gz
perl generate_training_input.pl mouse.median_expr.txt | gzip -c >Mouse_FantomAnnotations.InputData.pM10Kb.txt.gz
perl generate_training_input.pl Ouyang_mouseESC_RPKM_ensemblID.txt | gzip -c >Mouse_FantomAnnotations.InputData.pM10Kb.mESC.txt.gz

python setup_training_files.py -t 1000 -v 1000 Roadmap_FantomAnnotations.InputData.pM10Kb.txt.gz pM10Kb_1KTest
python setup_training_files.py -t 1000 -v 1000 Mouse_FantomAnnotations.InputData.pM10Kb.txt.gz pM10Kb_1KTest_Mouse

#run hyperparameter search for ~1000 iterations on a fast GPU (can take ~1-2 days to run)
python Xpresso.py tune tpe_1K_10epochs_optimized_0to20K.hyperopt pM10Kb_1KTest/

Rscript Fig1B_S2A.R FILE1 FILE2

for x in {1..10}; do { python Xpresso.py test tpe_1K_10epochs_optimized_0to20K.hyperopt pM10Kb_1KTest/ >trial_$x.txt; } done &

Rscript Fig1C.R FILE1
