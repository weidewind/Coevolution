use Math::Random::ISAAC;
#Module provides funktion for manipulation with observation vectors
#Observation vector is a vector of references to arrays of 3 elements (STR,X,Y):
#		STR - string is a text label, e.g. name of a tree branch
#		X - float
#		Y - uint is a number of observation of corresponding X value. Y=Y(X) is assumed
#	Observation vector has to be sorted by X

sub make_observation_vector{
	my ($rh_x,$rh_y)=@_;
	my @obsv;
	foreach my $name(keys %{$rh_x}){
		if(defined $rh_y->{$name}){
			next unless $rh_x->{$name}=~/^[\d\.+\-eE]+$/;
			next unless $rh_y->{$name}=~/^\d+$/;
			push @obsv, [($name,$rh_x->{$name},$rh_y->{$name})];
		};
	}
	@obsv=sort {$a->[1]<=>$b->[1]} @obsv;
	return @obsv;
}

#Function returns samples from CDF corresponding to the observation vector
#Returns hash where keys are indexes in the given observation vector, values are number of samples (<n)
sub sample_from_obsv{
	my $ra_obsv=shift;
	my $nsamples=shift;
	my $n=0;
	my %sample_idx;
	my @CPF;
	push @CPF,[@{$ra_obsv->[0]}[1,2]];
	for(my $i=1;$i<@{$ra_obsv};$i++){
		push @CPF,[@{$ra_obsv->[$i]}[1,2]];
		$CPF[-1]->[1]+=$CPF[-2]->[1];
	};
	die "\nError in sample_from_obsv(): Unable to throw $nsamples different samples fromthe total $CPF[-1]->[1] observations!" unless $nsamples<$CPF[-1]->[1];
	while($n<$nsamples){
		#my $rnum= int(rand($CPF[-1]->[1]));
		my $rnum= $rng->irand()%$CPF[-1]->[1];
		my $i=0;
		for(;$i<@CPF;$i++){
			last if($CPF[$i]->[1]>=$rnum)
		};
		if($ra_obsv->[$i]->[2]>$sample_idx{$i}){
			$sample_idx{$i}++;
			$n++;
		}
	}
	return %sample_idx;
};

#Generates a new randon observation vector with the same number of observations as in given vector
#	that are distributed proportionally to X values
sub shuffle_obsv{
	my $ra_obsv=shift;
	my @new_obsv;
	my $N=$ra_obsv->[0]->[2];
	$new_obsv[0]=[@{$ra_obsv->[0]}];
	$new_obsv[0]->[2]=0;
	for(my $i=1;$i<@{$ra_obsv};$i++){
		$N+=$ra_obsv->[$i]->[2]; # total number of observations
		$new_obsv[$i]=[@{$ra_obsv->[$i]}]; # copy $ra_obsv->[$i]
		$new_obsv[$i]->[2]=0;
		$new_obsv[$i]->[1]+=$new_obsv[$i-1]->[1]; # building "CDF" ("distribution") (sum not equal to one)
	};
	my %restrictor;
	while($N){
		#my $val=rand()*$new_obsv[-1]->[1];
		my $val=$rng->rand()*$new_obsv[-1]->[1]; # random value from 0 to sum of all x values
		my $i=0;
		for(;$i<@{$ra_obsv};$i++){					# select the least i at which "CDF" reaches $val (
			last if ($val<=$new_obsv[$i]->[1]); 
		};
		last if $restrictor{$i};					# prohibit addition of more than one observation
#die "\nOut of range" if $i==@{$ra_obsv};
		$new_obsv[$i]->[2]++;						# add observation to the corresponding array ( cannot be more than one in our case)
		$restrictor{$i} = 1;											
		$N--;
	};
	for(my $i=0;$i<@{$ra_obsv};$i++){
		$new_obsv[$i]->[1]=$ra_obsv->[$i]->[1];
#print "\n$new_obsv[$i]->[0]\t$new_obsv[$i]->[1]\t$new_obsv[$i]->[2]\t$ra_obsv->[$i]->[2]";
	};
	return @new_obsv;
}

sub mean_for_obsv{
	my $ra_obsv=shift;
	my $n=0;
	my $M=0;
	for(my $i=0;$i<@{$ra_obsv};$i++){
		$M+=$ra_obsv->[$i]->[1]*$ra_obsv->[$i]->[2];
		$n+=$ra_obsv->[$i]->[2];
	}
	return -1 if($n==0);
	return $M/$n;
};

sub median_for_obsv{
	my $ra_obsv=shift;
	my $n=0;
	my $N=0;
	for(my $i=0;$i<@{$ra_obsv};$i++){
		$N+=$ra_obsv->[$i]->[2];
	};
	my $i=0;
	for(;$i<@{$ra_obsv};$i++){
		$n+=$ra_obsv->[$i]->[2];
		last if $n>=$N/2;
	}
	return $ra_obsv->[$i]->[1];
};
