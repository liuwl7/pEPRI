#!/usr/bin/perl
my $sam_input=shift;
my $enhancer_bed=shift;
my $promoter_bed=shift;


my @interaction_files=(
$sam_input
);

foreach my $sam (@interaction_files){
	my $prefix=$sam;
	$prefix=~s/.+\///;
	$prefix=~s/.Chimeric.sam//;
	
	`perl from_sam_to_pair_reads_bed.pl $sam`;
	`bedtools intersect -wa -wb -a $enhancer_bed -b read_1.bed -F 0.5 > enhancer_overlap_with_read1.bed`;
	`bedtools intersect -wa -wb -a $enhancer_bed -b read_2.bed -F 0.5 > enhancer_overlap_with_read2.bed`;
	`bedtools intersect -wa -wb -a $promoter_bed -b read_1.bed -F 0.5 > promoter_overlap_with_read1.bed`;
	`bedtools intersect -wa -wb -a $promoter_bed -b read_2.bed -F 0.5 > promoter_overlap_with_read2.bed`;
	`perl collect_enhancer_promoter.pl enhancer_overlap_with_read1.bed enhancer_overlap_with_read2.bed promoter_overlap_with_read1.bed promoter_overlap_with_read2.bed $prefix`;
}