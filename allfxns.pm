package allfxns;

use Getopt::Long;
use POSIX qw/ceil floor/;
use List::Util qw/min max/;
#use Math::CDF qw/qnorm/;
#use List::MoreUtils qw/uniq/;
use Env; Env::import();

@ISA = qw(Exporter);
@EXPORT = qw(ceil floor qnorm uniq min max
		 bsub qsub parallelize 
		 fisher_yates_shuffle histogram log2 log10 mean median medianabsdev quantile round stdev sum trimmed_mean trimmed_stdev vecsum zscore 
		 alifold com comRNA fold gen_all_nmers readFasta rev revCom revComRNA gcContent cpgContent 
		 intersect unique array_diff array_minus);

#JOB SUBMISSION SUBROUTINES

sub bsub{
	($options, $error, $output, $job) = @_;
	system "bsub $options -e $error -o $output <$job";
	unlink $job;
}

sub qsub{
	($options, $error, $output, $job) = @_;
	system "qsub $options -S /bin/bash -e $error -o $output $job";
	unlink $job;
}

sub parallelize{
	%o = %{$_[0]};
	($cmd, $indir, $uniq, $instr, $outstr, $count, $jobnum) = ($o{"cmd"}, $o{"indir"}, $o{"uniq"}, $o{"instr"}, $o{"outstr"}, 0, 1000000);
	$outdir = $o{"outdir"} || "$TMP";
	$skipuniq = $o{"skipuniq"} || "";
	$exten = $o{"exten"} || "out";
	$qsub = $o{"qsub"} || 0;
	$bsub = $o{"bsub"} || 0;
	$onefile = $o{"onefile"} || 0;
	$skipsame = $o{"skipsame"} || 0;
	$redir = ($onefile) ? ">>" : ">";
	@files =  <$indir/*.$uniq>;
	system "mkdir $outdir" if !(-d $outdir);

	if ($qsub || $bsub) {
		$cmdperjob = $o{"cmdperjob"} || 1;
		$subopts = $o{"subopts"} || "-q idle";
		$suberr = $o{"suberr"} || "$TMP/sub.err";
		$subout = $o{"subout"} || "$TMP/sub.out";
		$jobfolder = $o{"jobfolder"} || "$TMP";
	}

	foreach $file (@files){
		++$count;
		$code = (split /\.$uniq/, (split /\//, $file)[-1])[0];
#		$fileexists = `grep -P '$code\\t' /lab/bartel3_ata/agarwal/metazoans/human.utrs/three_prime_UTR/bins.txt`; next if $fileexists;
#		$region = (split /\//, $file)[-2];
		$code = "concat$jobnum" if $onefile;
		next if (($skipsame && -s "$outdir/$code.$exten" != 0) || ($skipuniq ne "" && $file =~ /$skipuniq/)); # -s if file has zero size
		if (!$qsub && !$bsub) {
			#print "$cmd $instr $file $outstr $outdir/$code.$exten\n";
			#system "$cmd $instr $file $outstr $outdir/$code.$exten";
			#$numfile = `grep '>' $file | wc -l`;
			$outfile = $code;
			$numfile1 = (-e "$outdir/$outfile.$exten") ? `wc -l $outdir/$outfile.aln` : 1;
			print "$cmd $outdir/$outfile.aln\n" if $numfile1 == 0;
		}
		else {
			if ($cmd =~ /bin_MSA/ && $cmd !~ /all/){
#				($bin) = (split /\s/, `grep -m1 -P '^$code\t' $DIR/targetpred/robin/3UTRs_nonredundant_18577genes.UTR_cons.10bins`)[-1]; die "getbin" if $bin != int($bin); # $bin++; --> do this if robin's
				($tmp, $species, $region, $kmerlen) = (split /\s/, $cmd);
				($bin) = (split /\s/, `grep -m1 -P '^$code\t' $DIR/metazoans/$species/$region/allgenes.bins`)[-1] if $kmerlen == 2 || $kmerlen == 8; #.23way
#				print "$code, $species, $region, $kmerlen, $bin\n";
			}
##			print "$cmd $bin $instr $file 2>&- $outstr $redir $outdir/$code.$exten\n";
			$jobfile = "$jobfolder/job$jobnum.sh";
			open SH, ">>$jobfile" or die "can't open $jobfile";
			print SH "$cmd $bin $instr $file 2>&- $outstr $redir $outdir/$code.$exten\n";
			close SH;
			if ($count % $cmdperjob == 0) {
				qsub($subopts, $suberr, $subout, $jobfile) if $qsub;
				bsub($subopts, $suberr, $subout, $jobfile) if $bsub;
				$jobnum++;
			}
		}
	}
	qsub($subopts, $suberr, $subout, $jobfile) if $qsub;
	bsub($subopts, $suberr, $subout, $jobfile) if $bsub;
	print STDERR "COMPLETE!\n";
}

#STATISTICS SUBROUTINES

sub fisher_yates_shuffle {
	local $x = shift;
	for ($i = @$x; --$i; ) {
		$j = int rand ($i+1);
		next if $i == $j;
		@$x[$i,$j] = @$x[$j,$i];
	}
}

sub histogram {
	$bin_width = 10;
	if ($#_ == 1){
		($hash, $bin_width) = @_;
	}
	else { $hash = shift; }
	$max, $min;
	%a = %$hash;
	%histogram;
	foreach (keys %a){
		$histogram{ceil(($_ + 1) / $bin_width) -1} += $a{$_};
	}
	
	while ( ($key, $value) = each(%histogram) ) {
		$max = $key if !defined($min) || $key > $max;
		$min = $key if !defined($min) || $key < $min;
	}
	
	for ($i = $min; $i <= $max; $i++) {
		$bin = sprintf("% 10d", ($i) * $bin_width);
		$frequency = $histogram{$i} || 0;

		$frequency = "#" x $frequency;
		print $bin." ".$frequency."\n" if $frequency ne "";
	}

	print "===============================\n\n";
	print "		Width: ".$bin_width."\n";
	print "		Range: ".$min."-".$max."\n\n";
}

sub log2 { return log($_[0])/log(2); }

sub log10 { return log($_[0])/log(10); }

sub mean  { return sum($_[0])/scalar(@{$_[0]}); }

sub median{ return quantile($_[0], 2); }

sub medianabsdev{
	local $med = median($_[0]); local @b;
	push(@b, abs($_ - $med)) for (@{$_[0]});
	return median(\@b);
}

sub quantile{
	local $rpole = shift;
	local $x = shift;
	@pole = @$rpole;
	$ret;
	@pole = sort {$a <=> $b} @pole;
	if( ($#pole % $x) == 0 ) {
		$ret = $pole[int($#pole/$x)];
	} else {
		$ret = ($pole[int($#pole/$x)] + $pole[int($#pole/$x)+1]) / 2;
	}
	return $ret;
}

sub round{ return int($_[0] + 0.5 * ($_[0] <=> 0)); }

sub stdev{
      return 0 unless @_ > 1;
      local $mean = mean(\@_);
	local $tot = 0;
	foreach (@_){ $tot += ($_ - $mean)**2; }
      return sqrt( $tot / $#_ );
}

sub trimmed_mean {
	local $a = shift; local $perc = shift;
	$perc /= 200;
	@a = sort {$a <=> $b} @$a; $num = scalar(@a);
	@a = @a[int($num*$perc)..int($num*(1-$perc))];
	return mean(\@a);
}

sub trimmed_stdev {
	local $a = shift; local $perc = shift;
	$perc /= 200;
	@a = sort {$a <=> $b} @$a; $num = scalar(@a);
	@a = @a[int($num*$perc)..int($num*(1-$perc))];
	return stdev(@a);
}

sub sum{
	local $sum;
	$sum += $_ for(@{$_[0]});
	return $sum;
}

sub vecsum{	#sum two vectors passed as array refs
	$len = max(scalar(@{$_[0]}), scalar(@{$_[1]}))-1;
	for $i (0..$len){ ${$_[0]}[$i] += ${$_[1]}[$i]; }
}

sub zscore{
	local $val = shift;
	return ( $val - mean(\@_) ) / stdev(@_);
}


#NUCLEIC ACID SEQUENCE SUBROUTINES

sub alifold{
	local $file = shift;
	local $o = "-d0 -r -cv 0.6 -nc 0.5";
	$score = `RNAalifold $o $file | tail -1`;
	($mfe) = ($score =~ /.*\(\s*(-\d+.\d+) = .*\).*/);
	return $mfe;
}

sub com{
	local $seq = shift;
	$seq =~ tr/tucgaTUCGA/aagctAAGCT/;
	return $seq;
}

sub comRNA{
	local $seq = shift;
	$seq =~ tr/tucgaTUCGA/aagcuAAGCU/;
	return $seq;
}

sub fold{
	local $seq = shift;
	local $o = shift;
	local $score = `echo $seq | RNAfold $o | tail -1`;
	($mfe) = ($score =~ /.*\((\s?-\d+.\d+)\).*/);
	return $mfe;
}

sub gen_all_nmers{
	local $nmer_size = shift;
	local @words = '';
	foreach (1..$nmer_size) {
		@new_words = ();
		foreach $word (@words){
			foreach $i ( qw/A C T G/ ){ push (@new_words, $word.$i); }
		}
		@words = @new_words;
	}
	return @words;
}

sub readFasta{
	local $fasta = shift;
	local %fasta = ();
	open DNA, "<$fasta" || die "Could not open fasta file for $fasta\n";
	while ($line = <DNA>){ chomp $line;
		if ($line =~ /^>\s?(\w+\.?\d*)/){ $header = $1; }
#		if ($line =~ /^>(.*)/){ $header = $1; }
		else { $fasta{$header} .= $line; }
	}
	close DNA;
	return \%fasta;
}

sub rev{ return scalar reverse $_[0]; }

sub revCom{ return rev(com($_[0])); }

sub revComRNA{ return rev(comRNA($_[0])); }

sub gcContent{ return 0 if (() = ($_[0] =~ /[AUTCG]/ig)) == 0; return sprintf("%.3f", (() = ($_[0] =~ /[CG]/ig))/(() = ($_[0] =~ /[AUTCG]/ig))); } # %G + C in seq, ignoring case, missing nucleotides, or gaps

sub cpgContent{ return 0 if (() = ($_[0] =~ /[AUTCG]/ig)) == 0; return sprintf("%.3f", (() = ($_[0] =~ /CG/ig))/(length($_[0])-1)); }

#ARRAY SUBROUTINES

sub unique(@) {
	return keys %{ {map { $_ => undef } @_}}; 
}

sub intersect(\@\@) {
	my %e = map { $_ => undef } @{$_[0]};
	return grep { exists( $e{$_} ) } @{$_[1]};
}

sub array_diff(\@\@) {
	my %e = map { $_ => undef } @{$_[1]};
	return @{[ ( grep { (exists $e{$_}) ? ( delete $e{$_} ) : ( 1 ) } @{ $_[0] } ), keys %e ] };
}

sub array_minus(\@\@) {
	my %e = map{ $_ => undef } @{$_[1]};
	return grep( ! exists( $e{$_} ), @{$_[0]} ); 
}
