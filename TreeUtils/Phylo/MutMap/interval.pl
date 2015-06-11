#!/usr/bin/perl
use Statistics::Descriptive;

#percentiles();

sub percentiles {
open INTS, "<h3_syn_shuffler_10000_1";
open INTS_OUT, ">h3_syn_ints_nobins";
print  INTS_OUT "distance,mean,95_percentile,05_percentile\n";
my $t;
while(<INTS>){
	if ($_ =~ m/^interval/ || $t =~ m/^interval/){
		if ( $t =~ m/^interval/) {
			$_ = $t;
		}
		my @array;
		chomp($_);
		$_ =~ m/^interval:\s(\d+)/;
		print INTS_OUT $1.",";
		$t = <INTS>;
		chomp($t);
		push @array, $t;
		#print INTS_OUT $t."\n";
		while ($t = <INTS>){
			if (!($t =~ /^[\d-]/)){last;}
			chomp($t);
			push @array, $t;
			#print INTS_OUT $t." first\n";
		}
	#	for (my $i = 0; $i<3; $i++){
	#		my $counter = 0;
	#		while ($t = <INTS>){
	#		if (!($t =~ /^[\d-]/)){last;}
	#		chomp($t);
	#		$array[$counter] += $t;
	#		$counter++;
	#		#print INTS_OUT $t." added$i\n";
	#	}
	#	}
		#print INTS_OUT " now stats! \n";
		$stat = Statistics::Descriptive::Full->new();
  		$stat->add_data(\@array);
  		 print INTS_OUT $stat->mean();
  		print INTS_OUT ",";
  		print INTS_OUT  $stat->percentile(95);
  		print INTS_OUT ",";	
  		print INTS_OUT  $stat->percentile(5);
  		print INTS_OUT "\n";
	}

}
close INTS;
close INTS_OUT;
}

#trans_intervals();

  sub trans_intervals{
open INTS, "<n2_syn_shuffler_10000_1";
open INTS_OUT, ">n2_syn_shuffler_10000_1_trans";
my $t;
while(<INTS>){
	if ($_ =~ m/^interval/ || $t =~ m/^interval/){
		if ( $t =~ m/^interval/) {
			$_ = $t;
		}
		
		$t = <INTS>;
		chomp($t);
		print INTS_OUT $t.",";
		while ($t = <INTS>){
			if (!($t =~ /^[\d-]/)){last;}
			chomp($t);
			print INTS_OUT $t.",";
	}
  		print INTS_OUT "\n";
	}
}
close INTS;
close INTS_OUT;
  }
  
  collect_stats();
  sub collect_stats{
  	open IN, "<n2_pval";
  	open OUT, ">n2_bin_stats";
  	my $t;
  	while(<IN>){
  		$t = $_;
  		START: if ($t =~ m/^site\snumber\s([\d]+)/){
  			print OUT $1."\t";
  			$t = <IN>;
  			if ($t =~ m/^interval:\s[\d]+\s+([-e\.\d]+)/){
  				print OUT $1."\t";
  				$t = <IN>;
  				$t =~ m/^interval:\s[\d]+\s+([-e\.\d]+)/;
  				print OUT $1."\n";
  			}
  			else {
  				print OUT "\n";
  				goto START;
  			}
  		}
  	}
  }
  