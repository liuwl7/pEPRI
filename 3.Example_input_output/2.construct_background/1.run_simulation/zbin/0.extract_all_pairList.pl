#!/usr/bin/perl
use List::Util qw(shuffle); 
die "perl $0 ../0.source/HeLa_merge.noBackground.RRI.multiple.details.CountMultiple.network" if(@ARGV != 1);
my $net=shift;

my %network;
open(NT,$net) || die;
while(my $line=<NT>){
	chomp $line;
	my @sub=split/\s+/,$line;
        my $gene_a=$sub[0]."\t".$sub[1]."\t".$sub[2]."\t".$sub[3]."\t".$sub[4]."\t".$sub[5];
        my $gene_b=$sub[6]."\t".$sub[7]."\t".$sub[8]."\t".$sub[9]."\t".$sub[10]."\t".$sub[11];
        my $num=$sub[12];
        my @pair=($gene_a,$gene_b);
        @pair=sort @pair;
	$network{$pair[0]."\t".$pair[1]}=$num;
}


#check
my @all_pairs=sort {$network{$b} <=> $network{$a}} keys %network;
open(AP,">all_pairs.list") || die;
foreach (@all_pairs){
	print AP $_,"\t",$network{$_},"\n";
}
close AP;
