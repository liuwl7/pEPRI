#!/usr/bin/perl
open(IN,"result1.significant_interaction.interMole.network")||die;
open(OUT,">HNRNPU_merge.significant_interaction.interMole_E-P.network")||die;
my $ee=0;
my $ep=0;
my $pp=0;
while (my $line=<IN>){
	chomp $line;
	my @array=split/\s+/,$line;
	if(($array[3] =~ /E/ and $array[9] =~ /P/) or ($array[3] =~ /P/ and $array[9] =~ /E/)){
		$ep++;
		print OUT $line,"\n";

	}elsif($array[3] =~ /E/ and $array[9] =~ /E/){
		$ee++;
	}elsif($array[3] =~ /P/ and $array[9] =~ /P/){
		$pp++;
	}else{
		die;
	}

}
#print "ee",$ee,"\n";
#print "ep",$ep,"\n";
#print "pp",$pp,"\n";
