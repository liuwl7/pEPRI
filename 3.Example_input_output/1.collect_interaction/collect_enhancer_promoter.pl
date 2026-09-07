#!/usr/bin/perl
die "perl $0 enhancer_overlap_with_read1.bed enhancer_overlap_with_read2.bed promoter_overlap_with_read1.bed promoter_overlap_with_read2.bed\n" if(@ARGV != 5);
my $enhancer_overlap_read1=shift;
my $enhancer_overlap_read2=shift;
my $promoter_overlap_read1=shift;
my $promoter_overlap_read2=shift;
my $outputprefix=shift;

my %all_read_id;
my %read1_enhancer=read_overlap_files($enhancer_overlap_read1,1);
my %read2_enhancer=read_overlap_files($enhancer_overlap_read2,1);
my %read1_promoter=read_overlap_files($promoter_overlap_read1,2);
my %read2_promoter=read_overlap_files($promoter_overlap_read2,2);

my %effective_reads;
my %enhancer_promoter_interaction;
foreach my $id (keys %all_read_id){
	
	my @promoter_have_read1=keys %{$read1_promoter{$id}};
	my @enhancer_have_read1=keys %{$read1_enhancer{$id}};
	my @promoter_have_read2=keys %{$read2_promoter{$id}};
	my @enhancer_have_read2=keys %{$read2_enhancer{$id}};

	if(@promoter_have_read1 > 1 or @enhancer_have_read1 > 1 or @promoter_have_read2 > 1 or @enhancer_have_read2 > 1){
		print $id;
		die;
	}

	#check
	my @e_or_p_have_read1=(@promoter_have_read1,@enhancer_have_read1);
	my @e_or_p_have_read2=(@promoter_have_read2,@enhancer_have_read2);
	foreach my $l (@e_or_p_have_read1){
		foreach my $r (@e_or_p_have_read2){
			my @pair=($l,$r);
			@pair=sort @pair;
			$enhancer_promoter_interaction{$pair[0]."\t".$pair[1]}++;
			$effective_reads{$pair[0]."\t".$pair[1]}{$id}=1;
		}
	}
}

open(READS,">result1.$outputprefix.enhancer_link_promoter.chimeric_reads_id.list") || die;
foreach my $ep_pair (keys %effective_reads){
	foreach my $id (keys %{$effective_reads{$ep_pair}}){
		print READS $id,"\t",$ep_pair,"\n";
	}
}

open(OUTB,">result2.$outputprefix.enhancer_promoter_interaction.pair") || die;
foreach my $ep_pair (sort {$enhancer_promoter_interaction{$b} <=> $enhancer_promoter_interaction{$a}} keys %enhancer_promoter_interaction){
	print OUTB $ep_pair,"\t",$enhancer_promoter_interaction{$ep_pair},"\n";
}

sub read_overlap_files{
	my $file=shift;
	my $type=shift;
	my %hash;
	open(IN,$file) || die;
	while(my $line=<IN>){
		chomp $line;
		my @sub=split/\s+/,$line;
		if($type==1){
			my $enhancer=join"\t",@sub[0..3];
			$hash{$sub[7]}{$enhancer}=1;
			$all_read_id{$sub[7]}=1;
		}
		elsif($type==2){
			my $promoter=join"\t",@sub[0..3];
			$hash{$sub[7]}{$promoter}=1;
			$all_read_id{$sub[7]}=1;
		}
		else{
			die "wrong file type: $type\n";
		}
	}
	return %hash;
}

