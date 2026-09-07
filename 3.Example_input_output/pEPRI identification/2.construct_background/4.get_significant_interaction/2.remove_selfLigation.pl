#!/usr/bin/perl
die "perl $0 result1.significant_interaction.network\n" if(@ARGV != 1);
my $network=shift;

my $prefix=$network;
$prefix=~s/.network$//;

open(SL,">$prefix.intraMole.network") || die;
open(INTER,">$prefix.interMole.network") || die;
open(NT,$network) || die;
while(my $line=<NT>){
	chomp $line;
	my @sub=split/\s+/,$line;
	if($sub[3] eq $sub[9]){
		print SL $line,"\n";
	}
	else{
		print INTER $line,"\n";
	}
}
