#!/usr/bin/perl
open(SH,">run.sh") || die;
foreach my $t (1..20){
	print SH "nohup perl split_and_transpose_large_matrix.pl ../../../1.run_simulation/HNRNPU_merge/1.bacth_0-10000/result1.thread$t.simulation.list 100 thread$t &\n"
}
