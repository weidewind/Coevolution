
my $pdb_file = "C:/Users/weidewind/Documents/CMD/Coevolution/Influenza/2YP3_h3n2_with_receptor_analogue.pdb";
my $distance_file = "C:/Users/weidewind/Documents/CMD/Coevolution/Influenza/2YP3_pocket_distances.txt";
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

#in 1Ru7 numeration
my @pocket = (190, 191, 192, 193, 194, 195, 196, 197, 198, 135, 136, 137, 138, 221, 222, 223, 224, 225, 226, 227, 228, 98, 153, 183);	
my @neva_h1_pocket = (153,155,134,136,183,190,194,98,226,228);
my @neva_h3_pocket = (222,225,226,135,136,137,98,153,190,194);
my @n2_pocket = (151, 406, 277);
my @n2_bigger_pocket = (118,292,371,152,224,119,227,276,277,406);
my @n1_pocket = (151, 278, 402);
min_distance(\@n1_pocket);
sub min_distance{
	my $pdb_file = "C:/Users/weidewind/Documents/CMD/Coevolution/Influenza/PDB/4B7R.pdb";
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
 #	print $coord1[0]."\t".$coord2[0]."\t".$coord1[1]."\t".$coord2[1]."\t".$coord1[2]."\t".$coord2[2]."\n";
 	my $dist = sqrt(($coord1[0] - $coord2[0])**2 + ($coord1[1] - $coord2[1])**2 + ($coord1[2] - $coord2[2])**2);
 #	print $dist."\n";
 	return $dist;
 	
 }