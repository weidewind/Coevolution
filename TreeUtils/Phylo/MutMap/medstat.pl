sub logic_global_median_statistics{
	my $tree = parse_tree("C:/Users/weidewind/Documents/CMD/Coevolution/Influenza/Kryazhimsky11/h3.l.r.newick");
	my %fasta = parse_fasta("C:/Users/weidewind/Documents/CMD/Coevolution/Influenza/Kryazhimsky11/h3.all.fa");
	my @mutmaps = codonmutmap($tree, \%fasta);
	my %subs_on_node = %{$mutmaps[0]};
	my %nodes_with_sub = %{$mutmaps[1]};
	my $step = 1;
	my $iterate = 100;
	my @bootstrap_median_diff;
	my @bins;
	my @shuffler_bins;
#my @arr = (154,158,234,239);
print "site\tsame_median\tdiff_median\tmedian_difference\tpvalue\n";
#foreach my $ind(@arr){
	for (my $ind = 1; $ind <566; $ind++){
		print $ind."\t";
		my %distr = find_all_distances_except_seq_and_sis_radius($tree, \@{$nodes_with_sub{$ind}}, \%subs_on_node, $ind, 1);
		my @site_bins = distr_to_stathist(\%distr, $step);
		foreach $interval(@{$site_bins[0]}){
			$bins[0]->[$interval] += $site_bins[0]->[$interval];
		}
		foreach $interval(@{$site_bins[1]}){
			$bins[1]->[$interval] += $site_bins[1]->[$interval];
		}
	}
	my $same_median = hist_median(\@{$bins[0]});
	my $diff_median = hist_median(\@{$bins[1]});
	my $obs_difference = $diff_median-$same_median;
	print $same_median."\t"; #this have to be the median of "same" statistics
	print $diff_median."\t"; #this have to be the median of "diff" statistics
	print $obs_difference."\t";
	
	my @bootstrap_median_diff;
	for (my $t = 0; $t < $iterate; $t++){
		for (my $ind = 1; $ind <566; $ind++){
			my %shuffled_distr = shuffler($tree, \@{$nodes_with_sub{$ind}}, \%subs_on_node, $ind, 1);
			my @shuffler_site_bins = distr_to_stathist(\%shuffled_distr, $step);
			foreach $interval(@{$shuffler_site_bins[0]}){
				$shuffler_bins[0]->[$interval] += $shuffler_site_bins[0]->[$interval];
			}
			foreach $interval(@{$shuffler_site_bins[1]}){
				$shuffler_bins[1]->[$interval] += $shuffler_site_bins[1]->[$interval];
			}
			
		}
		push @bootstrap_median_diff, hist_median(\@{$shuffler_bins[1]})-hist_median(\@{$shuffler_bins[0]});
	}
	
	my @sorted_bootstrap = sort {$a <=> $b} @bootstrap_median_diff;
	my $pvalue = 0;
	for (my $i = 0; $i < $iterate; $i++){
		if($sorted_bootstrap[$i] >= $obs_difference){
			$pvalue = ($iterate - $i)/$iterate;
			last;
		}
	}
	print $pvalue."\n";

	
}