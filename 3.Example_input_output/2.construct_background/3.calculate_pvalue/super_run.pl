#!/usr/bin/perl
my @samples=(
"HNRNPU_merge",
);

foreach my $s (@samples){
	`mkdir $s`;
	`mkdir ./$s/1.bacth_0-10000/`;
	`mkdir ./$s/2.base_on_random/`;
	`cp ./zbin/*pl ./$s/1.bacth_0-10000/`;
	`cp ./zbin/*pl ./$s/2.base_on_random/`;
	open(WPL,">./$s/1.bacth_0-10000/write.pl")|| die;
	print WPL "#!/usr/bin/perl\n";
	print WPL "open(SH,\">run.sh\") || die;\n";
	print WPL "foreach my \$t (1..20){\n\tprint SH \"nohup perl pvalue_calculater.pl ../../../1.run_simulation/$s/1.bacth_0-10000/all_pairs.thread\$t.list ../../../2.pre-process/$s//1.bacth_0-10000/thread\$t.transposed.matrix thread\$t &\\n\"\n}\n";
	print WPL "print SH \"perl merge_pvalue.pl \";\n";
	print WPL "foreach my \$t (1..20){\n\tprint SH \"comparison.observed_to_simulated.thread\$t.xls \"\n};\n";
	print WPL "print SH \"> comparison.observed_to_simulated.finalMerge.xls\\n\";\n";
	close WPL;
	
	open(RWPL,">./$s/2.base_on_random/write.pl")|| die;
        print RWPL "#!/usr/bin/perl\n";
        print RWPL "open(SH,\">run.sh\") || die;\n";
        print RWPL "foreach my \$t (1..20){\n\tprint SH \"nohup perl pvalue_calculater.pl ../../../1.run_simulation/$s/2.base_on_random/all_pairs.thread\$t.list ../../../2.pre-process/$s/2.base_on_random/thread\$t.transposed.matrix thread\$t &\\n\"\n}\n";
        print RWPL "print SH \"perl merge_pvalue.pl \";\n";
        print RWPL "foreach my \$t (1..20){\n\tprint SH \"comparison.observed_to_simulated.thread\$t.xls \"\n};\n";
        print RWPL "print SH \"> comparison.observed_to_simulated.finalMerge.xls\\n\";\n";
        close RWPL;


	

}
