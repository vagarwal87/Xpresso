#region.1Mb.bed is human locus; region.1Mb.mouse.bed is mouse locus

BASEFILE="region.1Mb.intervals.100ntStep" #for Mouse locus use "region.1Mb.intervals.100ntStep.mouse"

#generate 100nt step with 10.5Kb window size
bedtools makewindows -b region.1Mb.bed -w 10500 -s 100 | perl -ne '@a=split/\t/; print $_ if $a[2]-$a[1] == 10500;' | sort | uniq >$BASEFILE.bed
#extract sequences from fasta of human or mouse genome
bedtools getfasta -tab -fi Homo_sapiens_assembly19.fasta -bed $BASEFILE.bed -fo $BASEFILE.input.txt

#generate predictions
python predict_seqs.py tpe_1K_10epochs_optimized_0to20K.hyperopt human_trainepoch.11-0.426.h5 $BASEFILE.input.txt $BASEFILE.Plus.txt
python predict_seqs.py --revCom tpe_1K_10epochs_optimized_0to20K.hyperopt human_trainepoch.11-0.426.h5 $BASEFILE.input.txt $BASEFILE.Minus.txt

#center plus and minus stranded predictions for given 100nt "TSS"
tail -n+2 $BASEFILE.Plus.txt | perl -ne '@a=($_ =~ /(.*):(\d+)-(\d+)\t(.*)/); print "chr".join("\t", $a[0], $a[1]+7000-5000, $a[1]+7000+5000, $a[3]+1), "\n";' >$BASEFILE.Plus.bedGraph
tail -n+2 $BASEFILE.Minus.txt | perl -ne '@a=($_ =~ /(.*):(\d+)-(\d+)\t(.*)/); print "chr".join("\t", $a[0], $a[1]+3500-5000, $a[1]+3500+5000, $a[3]+1), "\n";' >$BASEFILE.Minus.bedGraph

#convert to bigwig to view on UCSC browser
bedGraphToBigWig $BASEFILE.Plus.bedGraph hg19.chrom.sizes $BASEFILE.Plus.bw
bedGraphToBigWig $BASEFILE.Minus.bedGraph hg19.chrom.sizes $BASEFILE.Minus.bw
