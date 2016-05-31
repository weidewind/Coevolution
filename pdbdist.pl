use File::Spec;
use Cwd qw(abs_path cwd getcwd);

sub cox_distances {
	
	my $dir = getcwd();
	my $pdb_file  = File::Spec->catfile($dir, "1occ.pdb");
	my $pairs_file  = File::Spec->catfile($dir, "1occ_pairs_3.txt");
	my $distance_file  = File::Spec->catfile($dir, "1occ_distances_3.txt");

	open PDB, "<$pdb_file";
	while($line = <PDB>){
 		if ($line =~ /^ATOM\s+[0-9]+\s+CA\s/){
 			my @splitter = split (/\s+/, $line);
 			$distance_hash{$splitter[5].$splitter[4]} = [$splitter[6], $splitter[7], $splitter[8]];
 			print $splitter[5].$splitter[4]."\t";
 			print $splitter[6]."\t".$splitter[7]."\t".$splitter[8]."\n";
 		}
	}
	close PDB;
	open PAIRS,  "<$pairs_file";
	open OUT, ">$distance_file";
 	while($line = <PAIRS>){
 		chop($line);
 		my @splitter = split (/_/, $line);
 		my $bg = $splitter[0];
 		my $fg = $splitter[1];
 		my $distAA = distance($distance_hash{$bg."A"}, $distance_hash{$fg."A"});
 		my $distNN = distance($distance_hash{$bg."N"}, $distance_hash{$fg."N"});
 		my $distAN = distance($distance_hash{$bg."A"}, $distance_hash{$fg."N"});
 		my $distNA = distance($distance_hash{$bg."N"}, $distance_hash{$fg."A"});
 		my $samemin = ($distAA, $distNN)[$distAA > $distNN];
 		my $intermin = ($distAN, $distNA)[$distAN > $distNA];
 		my $min = ($intermin, $samemin)[$intermin > $samemin];
 		if ($samemin < $intermin){print OUT "same\t";}
 		else { print OUT "inter\t";}
 		print OUT $min."\n";
 	}
 	close PAIRS;
 	close OUT;
}

 sub distance(){	
 	my @coord1 = @{$_[0]};
 	my @coord2 = @{$_[1]};
 	my $dist = sqrt(($coord1[0] - $coord2[0])**2 + ($coord1[1] - $coord2[1])**2 + ($coord1[2] - $coord2[2])**2);
 	return $dist;
 }
 
cox_distances();