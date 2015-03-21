
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

 sub distance(){
 	
 	my @coord1 = @{$_[0]};
 	my @coord2 = @{$_[1]};
 #	print $coord1[0]."\t".$coord2[0]."\t".$coord1[1]."\t".$coord2[1]."\t".$coord1[2]."\t".$coord2[2]."\n";
 	my $dist = sqrt(($coord1[0] - $coord2[0])**2 + ($coord1[1] - $coord2[1])**2 + ($coord1[2] - $coord2[2])**2);
 #	print $dist."\n";
 	return $dist;
 	
 }