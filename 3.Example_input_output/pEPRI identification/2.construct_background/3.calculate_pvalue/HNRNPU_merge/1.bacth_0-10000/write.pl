#!/usr/bin/perl
open(SH,">run.sh") || die;
foreach my $t (1..20){
	print SH "nohup perl pvalue_calculater.pl ../../../1.run_simulation/HNRNPU_merge/1.bacth_0-10000/all_pairs.thread$t.list ../../../2.pre-process/HNRNPU_merge//1.bacth_0-10000/thread$t.transposed.matrix thread$t &\n"
}
print SH "perl merge_pvalue.pl ";
foreach my $t (1..20){
	print SH "comparison.observed_to_simulated.thread$t.xls "
};
print SH "> comparison.observed_to_simulated.finalMerge.xls\n";
