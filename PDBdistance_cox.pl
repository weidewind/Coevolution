


sub min_distance{
my $pdb_file = $_[0];
my $distance_file = "C:/Users/weidewind/Documents/CMD/Coevolution/Influenza/4B7R_pocket_distances.txt";
open PDB, "<$pdb_file" or die "$!";

my %distance_hash;
my %d_hash;

while($line = <PDB>){
 	if ($line =~ /^ATOM\s+[0-9]+\s+CA\s/){
 		my @splitter = split (/\s+/, $line);
 		$distance_hash{$splitter[5].$splitter[4]} = [$splitter[6], $splitter[7], $splitter[8]];
 	#	print $splitter[6]."\t".$splitter[7]."\t".$splitter[8]."\n";
 		
 	}
			}
close PDB; 
	open DIST, ">$distance_file";	
	my @pocket = @{$_[0]};		
for (keys %distance_hash){
	my $key1 = $_;

	my $num1 = $key1;
	chop($num1);
	for (@pocket){
		my $dist = distance($distance_hash{$key1}, $distance_hash{$_."A"});
		#my $num2 = $_;

		if (exists $d_hash{$num1}){
			if ($dist < $d_hash{$num1}){
				$d_hash{$num1} = $dist;
			}
		}
		else {
			$d_hash{$num1} = $dist;
		}
	}
}


for (keys %d_hash){
	print DIST $_."\t".$d_hash{$_}."\n";
}
close DIST;
}


sub pairs_distance {
	
}


sub cox_distances {
	
	my $pdb_file = "C:/Users/weidewind/workspace/perlCoevolution/cox/1occ.pdb";
	my $pairs_file = "C:/Users/weidewind/workspace/perlCoevolution/cox/1occ_pairs_3.txt";
	my $distance_file = "C:/Users/weidewind/workspace/perlCoevolution/cox/1occ_distances_3.txt";
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

cox_distances();

sub inter_distances {
	my $pdb_file = "C:/Users/weidewind/Documents/CMD/Coevolution/Influenza/1RU7_mod.pdb";
my $distance_file = "C:/Users/weidewind/Documents/CMD/Coevolution/Influenza/1RU7_mod_distances.txt";
open PDB, "<$pdb_file";

my %distance_hash;
my %d_hash;

while($line = <PDB>){
 	if ($line =~ /^ATOM\s+[0-9]+\s+CA\s/){
 		my @splitter = split (/\s+/, $line);
 		$distance_hash{$splitter[5].$splitter[4]} = [$splitter[6], $splitter[7], $splitter[8]];
 	#	print $splitter[6]."\t".$splitter[7]."\t".$splitter[8]."\n";
 		
 	}
			}
close PDB; 
open DIST, ">$distance_file";			
for (keys %distance_hash){
	my $key1 = $_;
	my $num1 = $key1;
	chop($num1);
	for (keys %distance_hash){
		my $dist = distance($distance_hash{$key1}, $distance_hash{$_});
		my $num2 = $_;
		chop($num2);
	#	print $num1."\t".$num2."\n";
		if (exists $d_hash{$num1."_".$num2}){
		#	print $key1."\t".$_."\n";
			if ($dist < $d_hash{$num1."_".$num2}){
			#	print $key1."\t".$_."\t".$dist."\n";
				$d_hash{$num1."_".$num2} = $dist;
			}
		}
		else {
			$d_hash{$num1."_".$num2} = $dist;
		}
	 #print DIST $_."\t".$key1."\t".$dist."\n";
	}
}


for (keys %d_hash){
	print DIST $_."\t".$d_hash{$_}."\n";
}
close DIST;

}
 sub distance(){	
 	my @coord1 = @{$_[0]};
 	my @coord2 = @{$_[1]};
 	my $dist = sqrt(($coord1[0] - $coord2[0])**2 + ($coord1[1] - $coord2[1])**2 + ($coord1[2] - $coord2[2])**2);
 	return $dist;
 }