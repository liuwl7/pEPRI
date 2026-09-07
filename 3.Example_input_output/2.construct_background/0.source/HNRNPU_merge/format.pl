#!/usr/bin/perl
die "perl $0 result4.enhancer_to_mergedPromoter.pair\n" if(@ARGV != 1);
my $enhancer_promoter_pair=shift;

open(EPP,$enhancer_promoter_pair) || die;
while(my $line=<EPP>){
	chomp $line;
	my @sub=split/\s+/,$line;
	foreach (0..3){
		print $sub[$_],"\t";
	}
	print "255\t+\t";
	foreach (4..7){
		print $sub[$_],"\t";
	}
	print "255\t+\t";
	print $sub[8],"\n";
}
	
