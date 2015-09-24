package DnaUtilities::compare;

use strict;
use vars qw(@ISA @EXPORT $VERSION);
use Exporter;
$VERSION = 1.00; # Or higher
@ISA = qw(Exporter);
@EXPORT = qw(count_substitutions nsyn_substitutions syn_substitutions nsyn_substitutions_codons is_neighbour_changing); # Symbols to autoexport (:DEFAULT tag)

use Bio::Tools::CodonTable;
use Class::Struct;

struct Substitution => {
	position => '$',
	ancestral_allele => '$',
	derived_allele => '$'
};

my $myCodonTable;
my %possible_subs;
my %neighbour_hash;

sub get_codon_table{
	if (!$myCodonTable){
		$myCodonTable = Bio::Tools::CodonTable->new();
	}
	return $myCodonTable
}

sub get_possible_subs {
	if (!%possible_subs){
		%possible_subs = (
			A => ['T', 'C', 'G'],
			T => ['A', 'C', 'G'],
			G => ['T', 'C', 'A'],
			C => ['T', 'A', 'G']
		);
	}
	return %possible_subs;
}


#This function counts synonimous and non synonymous substitutions between two sequences
sub count_substitutions{
	my $anc_seq=shift;
	my $seq=shift;
	my $rCDS=shift;
	my $ra_syn=shift;
	my $ra_nsyn=shift;
	@{$ra_syn}=();
	@{$ra_nsyn}=();
	return -1 unless defined($anc_seq)&&defined($seq);
	my $len=length $anc_seq;
	if(defined($rCDS)){
		$rCDS->[1]=$len-($len%3)-1 if $rCDS->[1]<0; 
	}else{
		$rCDS=[(0,$len-($len%3)-1)];
	};
	die "\nError in count_substitutions(): sequences have different lengths:\n$anc_seq\n$seq" if $len!=length $seq;
	my $myCodonTable   = Bio::Tools::CodonTable->new();
	for(my $i=$rCDS->[0];$i<=$rCDS->[1]-2;$i+=3){
		my $acod=substr $anc_seq,$i,3;
		my $cod=substr $seq,$i,3;
		next if $acod=~m/-/ || $cod=~m/-/;
		next if $acod eq $cod;
		next if $myCodonTable->is_ter_codon($acod) || $myCodonTable->is_ter_codon($cod);
		my $aaa=$myCodonTable->translate($acod);
		my $aa=$myCodonTable->translate($cod);
		$cod=substr $seq,1,$i;
		my $n=0;
		while($cod=~m/-/g){$n++};
		$n=($i-$rCDS->[0]-$n)/3+1;
		if($aaa eq $aa){
			push @{$ra_syn},$n;
		}else{
			my $p=Substitution->new();
			$p->position($n);
			$p->ancestral_allele($aaa);
			$p->derived_allele($aa);
			push @{$ra_nsyn},$p;
		};
	};
	return @{$ra_syn}+@{$ra_nsyn}; # sum of lengths
};



sub nsyn_substitutions{
	my $anc_seq=shift;
	my $seq=shift;
	my $rCDS=shift;
	#my $ra_nsyn=shift;
	my %ra_nsyn;
	#@{$ra_nsyn}=();
	#%{$ra_nsyn} = {};
	return -1 unless defined($anc_seq)&&defined($seq);
	my $len=length $anc_seq;
	if(defined($rCDS)){
		$rCDS->[1]=$len-($len%3)-1 if $rCDS->[1]<0; 
	}else{
		$rCDS=[(0,$len-($len%3)-1)];
	};
	die "\nError in count_substitutions(): sequences have different lengths:\n$anc_seq\n$seq" if $len!=length $seq;
	my $myCodonTable   = Bio::Tools::CodonTable->new();
	for(my $i=$rCDS->[0];$i<=$rCDS->[1]-2;$i+=3){
		my $acod=substr $anc_seq,$i,3;
		my $cod=substr $seq,$i,3;
		next if $acod=~m/-/ || $cod=~m/-/;
		next if $acod eq $cod;
		next if $myCodonTable->is_ter_codon($acod) || $myCodonTable->is_ter_codon($cod);
		my $aaa=$myCodonTable->translate($acod);
		my $aa=$myCodonTable->translate($cod);
		$cod=substr $seq,1,$i;
		my $n=0;
		while($cod=~m/-/g){$n++};
		$n=($i-$rCDS->[0]-$n)/3+1;
		if($aaa ne $aa){
			my $p=Substitution->new();
			$p->position($n);
			$p->ancestral_allele($aaa);
			$p->derived_allele($aa);
			#push @{$ra_nsyn},$p;
			print ("mutmap: ".$n."\t".$aa."\t".$aaa."\n");
			$ra_nsyn{$n} = $p;
		};
	};
	##return @{$ra_nsyn}; 
	return  %ra_nsyn;
}

sub nsyn_substitutions_codons{
	my $anc_seq=shift;
	my $seq=shift;
	my $rCDS=shift;
	#my $ra_nsyn=shift;
	my %ra_nsyn;
	#@{$ra_nsyn}=();
	#%{$ra_nsyn} = {};
	return -1 unless defined($anc_seq)&&defined($seq);
	my $len=length $anc_seq;
	if(defined($rCDS)){
		$rCDS->[1]=$len-($len%3)-1 if $rCDS->[1]<0; 
	}else{
		$rCDS=[(0,$len-($len%3)-1)];
	};
	die "\nError in count_substitutions(): sequences have different lengths:\n$anc_seq\n$seq" if $len!=length $seq;
	my $myCodonTable   = Bio::Tools::CodonTable->new();
	for(my $i=$rCDS->[0];$i<=$rCDS->[1]-2;$i+=3){
		my $acod=substr $anc_seq,$i,3;
		my $cod=substr $seq,$i,3;
		next if $acod=~m/-/ || $cod=~m/-/;
		next if $acod eq $cod;
		next if $myCodonTable->is_ter_codon($acod) || $myCodonTable->is_ter_codon($cod);
		my $aaa=$myCodonTable->translate($acod);
		my $aa=$myCodonTable->translate($cod);
		my $codprev=substr $seq,1,$i;
		my $n=0;
		while($codprev=~m/-/g){$n++};
		$n=($i-$rCDS->[0]-$n)/3+1;
		if($aaa ne $aa){
			my $p=Substitution->new();
			$p->position($n);
			$p->ancestral_allele($acod);
			$p->derived_allele($cod);
			#push @{$ra_nsyn},$p;
			#print ("mutmap: ".$n."\t".$aa."\t".$aaa."\n");
			$ra_nsyn{$n} = $p;
		};
	};
	##return @{$ra_nsyn}; 
	return  %ra_nsyn;
}

sub syn_substitutions{
	my $anc_seq=shift;
	my $seq=shift;
	my $rCDS=shift;
	#my $ra_nsyn=shift;
	my %ra_nsyn;
	#@{$ra_nsyn}=();
	#%{$ra_nsyn} = {};
	return -1 unless defined($anc_seq)&&defined($seq);
	my $len=length $anc_seq;
	if(defined($rCDS)){
		$rCDS->[1]=$len-($len%3)-1 if $rCDS->[1]<0; 
	}else{
		$rCDS=[(0,$len-($len%3)-1)];
	};
	die "\nError in count_substitutions(): sequences have different lengths:\n$anc_seq\n$seq" if $len!=length $seq;
	my $myCodonTable   = Bio::Tools::CodonTable->new();
	for(my $i=$rCDS->[0];$i<=$rCDS->[1]-2;$i+=3){
		my $acod=substr $anc_seq,$i,3;
		my $cod=substr $seq,$i,3;
		next if $acod=~m/-/ || $cod=~m/-/;
		next if $acod eq $cod;
		next if $myCodonTable->is_ter_codon($acod) || $myCodonTable->is_ter_codon($cod);
		my $aaa=$myCodonTable->translate($acod);
		my $aa=$myCodonTable->translate($cod);
		my $codprev=substr $seq,1,$i;
		my $n=0;
		while($codprev=~m/-/g){$n++};
		$n=($i-$rCDS->[0]-$n)/3+1;
		if($aaa eq $aa){
			my $p=Substitution->new();
			$p->position($n);
			$p->ancestral_allele($acod);
			$p->derived_allele($cod);
			#push @{$ra_nsyn},$p;
			#print ("mutmap: ".$n."\t".$acod."\t".$cod."\n");
			$ra_nsyn{$n} = $p;
		};
	};
	##return @{$ra_nsyn}; 
	return  %ra_nsyn;
}

# tells if synonimous substitution changes the range of one-symbol neighbours
sub is_neighbour_changing {
	my $subst = $_[0];
	my $full = $_[1];
	if (!$full){
		$full = 0;
	}
	
	my $existing_answer = $neighbour_hash{$subst->{"Substitution::ancestral_allele"}}->{$subst->{"Substitution::derived_allele"}}->{$full};
	if ($existing_answer){
		return $existing_answer
	}
	else {
	
	my $codonTable = get_codon_table();
	
	my %anc_neighbours = get_neighbours($subst->{"Substitution::ancestral_allele"});
	my %der_neighbours = get_neighbours($subst->{"Substitution::derived_allele"});
	
	my $answer = 0;
	if (length(keys %anc_neighbours) != length(keys %der_neighbours)){
		$answer = 1;
	}
	else {
		if ($full == 1){
			foreach my $k (keys %anc_neighbours){
				if (!$der_neighbours{$k} || $der_neighbours{$k} ne $anc_neighbours{$k}){
					$answer = 1;
				}
			}
		}
		else {
			foreach my $k (keys %anc_neighbours){
				if (!$der_neighbours{$k}){
					$answer = 1;
				}
			}
		}
	}
	
	$neighbour_hash{$subst->{"Substitution::ancestral_allele"}}->{$subst->{"Substitution::derived_allele"}}->{$full} = $answer;
	return $answer;
	}

}
test_is_neighbour_changing();

sub test_is_neighbour_changing{
			my $p=Substitution->new();
			$p->position(5);
			$p->ancestral_allele("CCT");
			$p->derived_allele("CCC");
			print is_neighbour_changing($p, 1);
			my $p=Substitution->new();
			$p->position(5);
			$p->ancestral_allele("CCA");
			$p->derived_allele("CCC");
			print is_neighbour_changing($p, 1);
			$p->position(5);
			$p->ancestral_allele("CCT");
			$p->derived_allele("CCC");
			print is_neighbour_changing($p, 1);
			$p->position(5);
			$p->ancestral_allele("gcg");
			$p->derived_allele("gcc");
			print is_neighbour_changing($p, 1);
}

sub get_neighbours {
	my $codon = $_[0];
	my %neighbours;
	my $codonTable = get_codon_table();
	my $possible_subs = get_possible_subs();
	for (my $i = 0; $i < 3; $i++){
		my $str = $codon;
		foreach my $letter (@{$possible_subs{substr($codon, $i, 1)}}){
			substr($str, $i, 1) = $letter;
			$neighbours{$codonTable->translate($str)}++;
		}
	}
	return %neighbours;
}

sub test_get_neighbours {
	my %neigh = get_neighbours("ATG");
	foreach my $n(keys %neigh){
		print $n."\t".$neigh{$n}."\n";
	}
}

#test_is_neighbour_changing();

1;