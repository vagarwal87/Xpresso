#!/usr/bin/perl

open IN, "<ensembl2geneName_v90_mm10.txt";
while(<IN>){ chomp;
	@a = split /\t/, $_;
	$id2parent{$a[1]} = $a[0] if $a[2] =~ /^\d+|^X/; #remove haplotypes when considering Ensembl ID
}
close IN;

open IN, "cut -f 2,10 MGI_EntrezGene.rpt | ";
while(<IN>){ chomp;
	@a = split /\t/, $_;
	$mgi2synonym{$a[0]} = $a[1];
}
close IN;

open IN, "cut -f 1,3 ensembl2entrezID_v90_mm10.txt | ";
while(<IN>){ chomp;
	@a = split /\t/, $_;
	$entrez2ensembl{$a[1]} = $a[0] if $a[1] ne '';
}
close IN;

open IN, "cut -f 3,11 MGI_Gene_Model_Coord.rpt | ";
while(<IN>){ chomp;
	@a = split /\t/, $_;
	$mgi2parent{$a[0]} = $a[1];
	@b = split /\|/, $mgi2synonym{$a[0]};
	foreach $i (@b){ $mgi2parent{$i} = $a[1]; }
}
close IN;

open IN, "cut -f 3,11 MGI_Gene_Model_Coord.rpt | ";
while(<IN>){ chomp;
	@a = split /\t/, $_;
	$mgi2parent{$a[0]} = $a[1];
}
close IN;

open IN, "<Ouyang_mESC_RPKM.txt";
open OUT, ">Ouyang_mESC_RPKM_ensemblID.txt";
while(<IN>){ chomp; @a=split; $id = $a[0];
	$id = $entrez2ensembl{$a[0]};
	$id = $id2parent{$a[1]} if $id eq '';
	print OUT join("\t", $id, $a[2]),"\n" if $id ne '' && !$seenid{$id};
	print STDERR $_,"\n" if $id eq '' || $seenid{$id};
	$seenid{$id} = 1;
}
close IN;
close OUT;

open IN, "<mm10_cage_promoters.bed";
open OUT, ">mm10_cage_promoters_ensemblID.bed";
while(<IN>){ chomp; @a=split; $id = $a[3];
	@ids = split /,/, $id;
	foreach $id (@ids){
		($promoter, $gene) = split /\@/, $id;
		$a[3] = $id2parent{$gene};
		$a[3] = $mgi2parent{$gene} if $a[3] eq '';
		print OUT join("\t", @a),"\n" if $promoter eq 'p1' && $a[3] =~ /ENS/ && !$seenid{$a[3]};
		$seenid{$a[3]} = 1;
	}
}
close IN;
close OUT;