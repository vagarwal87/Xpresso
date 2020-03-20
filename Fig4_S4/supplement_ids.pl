#!/usr/bin/perl

open IN, "<hg19_promoters_cage_corrected_withChr.bed";
while(<IN>){
	($id) = ($_ =~ /(ENSG\d+)/);
	$seen{$id}=1;
	print $_;
}
close IN;

open IN, "zgrep -P '\tgene\t' gencode.v27lift37.basic.annotation.gtf.gz | grep protein_coding | ";
while(<IN>){ chomp;
	@a = split /\t/;
	($id) = ($a[-1] =~ /(ENSG\d+)/);
	print join("\t", $a[0], $a[4]-1, $a[4]+1, $id, '0', $a[6]), "\n" if !$seen{$id} && $a[6] eq '-';
	print join("\t", $a[0], $a[3]-1, $a[3]+1, $id, '0', $a[6]), "\n" if !$seen{$id} && $a[6] eq '+';
	$seen{$id}=1;
}
close IN;