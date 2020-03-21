#!/usr/bin/perl

open IN, "<ensembl2geneName_v90.txt";
while(<IN>){ chomp;
	@a = split /\t/, $_;
	$id2parent{$a[1]} = $a[0] if $a[2] =~ /^\d+|^X/; #remove haplotypes when considering Ensembl ID
}
close IN;

open IN, "<HGNC2Ensembl.txt";
while(<IN>){ chomp;
	@a = split /\t/, $_;
	@b = split /, /, $a[1];
	@c = split /, /, $a[2];
	$hgnc2parent{$a[0]} = $a[3];
	foreach $i (@b){ $hgnc2parent{$i} = $a[3]; }
	foreach $i (@c){ $hgnc2parent{$i} = $a[3]; }
}
close IN;

open IN, "<hg38_cage_promoters.bed";
open OUT, ">hg38_cage_promoters_ensemblID.bed";
while(<IN>){ chomp; @a=split; $id = $a[3];
	@ids = split /,/, $id;
	foreach $id (@ids){
		($promoter, $gene) = split /\@/, $id;
		$a[3] = $id2parent{$gene};
		$a[3] = $hgnc2parent{$gene} if $a[3] eq '';
		print OUT join("\t", @a),"\n" if $promoter eq 'p1' && $a[3] ne '';
	}
}
close IN;
close OUT;