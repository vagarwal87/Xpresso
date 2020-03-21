#!/usr/bin/perl

$file = shift;
$spec = shift;

$dist = 10000;

open IN, "zgrep -P '\texon\t' $file | ";
while(<IN>){
	($chr, $region, $start, $stop, $str, $last) = (split /\t/)[0,2,3,4,6,-1];
	($parent) = ($last =~ /(ENS\w*G\d+)/);
	$has5pUTR{$id} = 1;
	if($str eq '+'){
		$allregions{$parent} = $start if (! exists $allregions{$parent} || $start < $allregions{$parent});
	}
	else{
		$allregions{$parent} = $stop if (! exists $allregions{$parent} || $stop > $allregions{$parent});
	}
}
close IN;

open IN, "zgrep -P '\texon\t' $file | ";
while(<IN>){
	($start, $stop, $str, $last) = (split /\t/)[3,4,6,-1];
	($parent) = ($last =~ /(ENS\w*G\d+)/);
	next if $seenids{$parent};
	next if $str eq '+' && $allregions{$parent} != $start;
	next if $str eq '-' && $allregions{$parent} != $stop;
	@a = split /\t/, $_;
	$a[-1] = $parent;
	$a[2] = $parent;
	$a[3] = $allregions{$parent} - $dist;
	$a[4] = $allregions{$parent} + $dist;
	if ($spec eq "mouse"){ print join("\t", 'chr'.$a[0], $a[3], $a[4], $a[-1], 0, $str), "\n" if $a[3] > 0; }
	else { print join("\t", $a[0], $a[3], $a[4], $a[-1], 0, $str), "\n" if $a[3] > 0; }
	$seenids{$parent} = 1;
}
close IN;