#!/usr/bin/perl
my @samples=(
"HNRNPU_merge",
);

foreach my $s (@samples){
	`mkdir $s`;
	`mkdir ./$s/1.bacth_0-10000/`;
	`mkdir ./$s/2.base_on_random/`;
	`cp ./zbin/1.MonteCarlo_simulation.pl ./$s/1.bacth_0-10000/`;
	`cp ./zbin/1.MonteCarlo_simulation.pl ./$s/2.base_on_random/`;
	`cp ./zbin/0.creat_random_interaction.pl ./$s/2.base_on_random/`;
	`cp ./zbin/0.extract_all_pairList.pl ./$s/2.base_on_random/`;
	open(WPL,">./$s/1.bacth_0-10000/write.pl")|| die;
	print WPL "open(SH,\">run.sh\") || die;\n";
	print WPL "foreach (1..20){\n\tprint SH \"nohup perl 1.MonteCarlo_simulation.pl ../../../0.source/$s/result3.$s.enhancer_to_promoter.format.pair 5000 thread\$_ > result1.thread\$_.simulation.list &\\n\"\n}";
	close WPL;
	
	open(RWSH,">./$s/2.base_on_random/work.sh") || die;
	print RWSH "perl 0.creat_random_interaction.pl ../../../0.source/$s/result3.$s.enhancer_to_promoter.format.pair > result0.random.network\n";
	close RWSH;
	open(RWPL,">./$s/2.base_on_random/write.pl")|| die;
	print RWPL "open(SH,\">run.sh\") || die;\n";
	print RWPL "foreach (1..20){\n\tprint SH \"nohup perl 1.MonteCarlo_simulation.pl result0.random.network 5000 thread\$_ > result1.thread\$_.simulation.list &\\n\"\n}";
	close RWPL;
	

}
