#!/usr/bin/perl
my @samples=(
"HNRNPU_merge",
);

foreach my $s (@samples){
	`mkdir $s`;
	`mkdir ./$s/1.bacth_0-10000/`;
	`mkdir ./$s/2.base_on_random/`;
	`cp ./zbin/split_and_transpose_large_matrix.pl ./$s/1.bacth_0-10000/`;
	`cp ./zbin/split_and_transpose_large_matrix.pl ./$s/2.base_on_random/`;
	open(WPL,">./$s/1.bacth_0-10000/write.pl")|| die;
	print WPL "#!/usr/bin/perl\n";
	print WPL "open(SH,\">run.sh\") || die;\n";
	print WPL "foreach my \$t (1..20){\n\tprint SH \"nohup perl split_and_transpose_large_matrix.pl ../../../1.run_simulation/$s/1.bacth_0-10000/result1.thread\$t.simulation.list 100 thread\$t &\\n\"\n}\n";
	close WPL;
	
	open(RWPL,">./$s/2.base_on_random/write.pl")|| die;
	print RWPL "#!/usr/bin/perl\n";
	print RWPL "open(SH,\">run.sh\") || die;\n";
	print RWPL "foreach my \$t (1..20){\n\tprint SH \"nohup perl split_and_transpose_large_matrix.pl ../../../1.run_simulation/$s/2.base_on_random/result1.thread\$t.simulation.list 100 thread\$t &\\n\"\n}\n";
	close RWPL;
	

}
