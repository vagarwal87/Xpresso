#!/usr/bin/perl

$species = shift;

if ($species eq "human"){ open IN, "zcat Roadmap_FantomAnnotations.InputData.pM10Kb.txt.gz | cut -f 1,11 | "; }
else{ open IN, "zcat Mouse_FantomAnnotations.InputData.pM10Kb.txt.gz | cut -f 1,11 | "; }
while(<IN>){
    ($id, $seq) = split /\t/, $_;
    $id2seq{$id} = $seq;
}
close IN;

while(<>){ chomp;
    @a = split /\t/, $_;
    $id = $a[0];
    $start = $a[1]+3000-1;
    $len = $a[3];
    print ">$id\_$start\_$len\n", substr($id2seq{$id}, $start, $len), "\n";
}
