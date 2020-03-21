#!/usr/bin/perl

use allfxns;

$exprMat = shift;
$spec = ($exprMat =~ /mouse/)? 1 : 0; #species is mouse or human?

sub readBed{
	local $bed = shift;
	local %bed = ();
	open BEDENTRY, "<$bed" || die "Could not open bed file for $bed\n";
	while ($line = <BEDENTRY>){
		@a = split /\t/, $line;
		$bed{$a[3]} = $line;
	}
	close BEDENTRY;
	return \%bed;
}

sub readFasta{
	local $fasta = shift;
	local %fasta = ();
	open DNA, "zcat $fasta | " || die "Could not open fasta file for $fasta\n";
	while ($line = <DNA>){ chomp $line;
		if ($line =~ /^>\s?(\w+\d)\.?\d*/){ $header = $1; }
		else { $fasta{$header} .= $line; }
	}
	close DNA;
	return \%fasta;
}

if ($spec){
	open IN, "zcat Mus_musculus.GRCm3.90.chosenTranscript.gtf.gz | ";
}
else{
	open IN, "zcat Homo_sapiens.GRCh38.90.chosenTranscript.gtf.gz | ";
}
while(<IN>){
	($region, $start, $stop, $last) = (split /\t/)[2,3,4,-1];
	($id) = ($last =~ /(ENS\w*GR?[\d|\.]+)/);
	$lengths{$id}{$region} += ($stop-$start);
	$cdsexoncount{$id}++ if $region eq 'CDS';
}
close IN;

if ($spec){
	%promoterbed = %{ readBed("mm10_promoters.bed") };
	%promoters = %{ readFasta("mm10_promoters.fa.gz") };
	%fantombed = %{ readBed("mm10_cage_promoters_ensemblID.bed") };
	%fantompromoters = %{ readFasta("mm10_cage_promoters_ensemblID.fa.gz") };
	%utr5p = %{ readFasta("mm10_ensembl90_5utrs.fa.gz") };
	%orfs = %{ readFasta("mm10_ensembl90_orfs.fa.gz") };
	%utr3p = %{ readFasta("mm10_ensembl90_3utrs.fa.gz") };
}
else{
	%promoterbed = %{ readBed("hg38_promoters.bed") };
	%fantombed = %{ readBed("hg38_cage_promoters_ensemblID.bed") };
	%promoters = %{ readFasta("hg38_promoters.fa.gz") };
	%fantompromoters = %{ readFasta("hg38_cage_promoters_ensemblID.fa.gz") };
	%utr5p = %{ readFasta("hg38_ensembl90_5utrs.fa.gz") };
	%orfs = %{ readFasta("hg38_ensembl90_orfs.fa.gz") };
	%utr3p = %{ readFasta("hg38_ensembl90_3utrs.fa.gz") };
}

if ($spec){
	open BED, ">mm10_promoters_cage_corrected.bed";
}
else{
	open BED, ">hg38_promoters_cage_corrected.bed";
}
print join("\t", "ENSID", "EXPRESSION", "UTR5LEN", "CDSLEN", "INTRONLEN", "UTR3LEN", "UTR5GC", "CDSGC", "UTR3GC", "ORFEXONDENSITY", "PROMOTER"), "\n";
open IN, "<$exprMat";
while(<IN>){ chomp;
	@a=split /\t/;
	$id = $a[0];
	if (($promoters{$id} ne '' || $fantompromoters{$id} ne '') && $lengths{$id}{"CDS"} > 0){
		$promoter = $promoters{$id};
		$promoterbed = $promoterbed{$id};
		if ($fantompromoters{$id} ne ''){
			$promoter = $fantompromoters{$id};
			$promoterbed = $fantombed{$id};
		}
		print BED $promoterbed;
		print join("\t", $id, $a[1], int($lengths{$id}{"five_prime_utr"}), int($lengths{$id}{"CDS"}), 
			int($lengths{$id}{"transcript"})-(int($lengths{$id}{"three_prime_utr"})+int($lengths{$id}{"CDS"})+int($lengths{$id}{"five_prime_utr"})), 
			int($lengths{$id}{"three_prime_utr"}), gcContent($utr5p{$id}), gcContent($orfs{$id}), gcContent($utr3p{$id}), 
			sprintf("%.2f", $cdsexoncount{$id}*1000/$lengths{$id}{"CDS"}), $promoter), "\n";
	}
	else{ $count++; }
}
close BED;

print STDERR "$count IDs missing/revised due to annotation version changes\n";