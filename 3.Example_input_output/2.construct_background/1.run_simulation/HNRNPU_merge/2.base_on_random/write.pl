open(SH,">run.sh") || die;
foreach (1..20){
	print SH "nohup perl 1.MonteCarlo_simulation.pl result0.random.network 5000 thread$_ > result1.thread$_.simulation.list &\n"
}