open(SH,">run.sh") || die;
foreach (1..20){
	print SH "nohup perl 1.MonteCarlo_simulation.pl ../../../0.source/HNRNPU_merge/result3.HNRNPU_merge.enhancer_to_promoter.format.pair 5000 thread$_ > result1.thread$_.simulation.list &\n"
}