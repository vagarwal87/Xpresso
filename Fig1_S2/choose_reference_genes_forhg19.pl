#!/usr/bin/perl

$file = shift;

open IN, "zgrep -P '\tUTR|CDS\t' $file | ";
while(<IN>){
	($region, $start, $stop, $last) = (split /\t/)[2,3,4,-1];
	($parent, $id) = ($last =~ /gene_id "(ENS\w*GR?[\d|\.]+.*)"; transcript_id "(ENS\w*TR?[\d|\.]+.*)"; gene_type/);
	$lengths{$id}{$region} += ($stop-$start);
}
close IN;

open IN, "zgrep -P '\tCDS\t' $file | ";
while(<IN>){
	($parent, $id) = ($_ =~ /gene_id "(ENS\w*GR?[\d|\.]+.*)"; transcript_id "(ENS\w*TR?[\d|\.]+.*)"; gene_type/);
	$reptranscript{$parent} = $id if (! defined $reptranscript{$parent});
	$reptranscript{$parent} = $id if ($lengths{$id}{"CDS"} > $lengths{$reptranscript{$parent}}{"CDS"} && $lengths{$id}{"UTR"} > 0);
}
close IN;

%okids = map { $_ => 1 } values %reptranscript;
foreach (keys %reptranscript){
	$repid{$reptranscript{$_}} = $_;
}

open IN, "zcat $file |";
while(<IN>){
	($chr, $region, $start, $stop, $str, $last) = (split /\t/)[0,2,3,4,6,-1];
	($parent, $id) = ($last =~ /gene_id "(ENS\w*GR?[\d|\.]+.*)"; transcript_id "(ENS\w*TR?[\d|\.]+.*)"; gene_type/);
	$repid = $repid{$id};
	if ($okids{$id}){
		@a = split /\t/, $_;
		$a[-1] = "$parent";
		$_ = join ("\t", @a)."\n";
		print "$_";
	}
}
close IN;