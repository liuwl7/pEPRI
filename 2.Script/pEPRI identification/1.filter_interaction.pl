#!/usr/bin/perl
die "perl $0 network contact_cutoff raw_pvalue_cutoff\n" if(@ARGV != 3);
my $net=shift;
my $contact_cutoff=shift;
my $raw_pvalue_cutoff=shift;

my %gene_targetNum;
my $significant_interactions;
my $nonSignificant_interactions;
open(NT,$net) || die;
while(my $line=<NT>){
	chomp $line;
	my @sub=split/\s+/,$line;
	my $num=$sub[12];
	my $raw_pvalue=$sub[14];
	my $gene_a=join"\t",@sub[0..5];
	my $gene_b=join"\t",@sub[6..11];

	if($num >= $contact_cutoff and $raw_pvalue <= $raw_pvalue_cutoff){
		$gene_targetNum{$gene_a}++;
		$gene_targetNum{$gene_b}++;
		$significant_interactions++;
		print $line,"\n";
	}
	else{
		$nonSignificant_interactions++;
	}
}

warn "$significant_interactions significant interactions\n";
warn "$nonSignificant_interactions other types interactions\n";
