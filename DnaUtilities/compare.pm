package DnaUtilities::compare;

use strict;
use vars qw(@ISA @EXPORT $VERSION);
use Exporter;
$VERSION = 1.00; # Or higher
@ISA = qw(Exporter);
@EXPORT = qw(count_substitutions, nsyn_substitutions); # Symbols to autoexport (:DEFAULT tag)

use Bio::Tools::CodonTable;
use Class::Struct;

struct Substitution => {
	position => '$',
	ancestral_allele => '$',
	derived_allele => '$'
};

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


1;