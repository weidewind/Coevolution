#!/usr/bin/perl
## This class provides a map of mutations (node_name -> array of sites that mutated on the branch),  
## given a tree and sequences of all nodes (including internal ones)

use strict;
use Bio::Phylo::IO;
use DnaUtilities::compare qw(nsyn_substitutions);
use Bio::SeqIO;
use Data::Dumper;
use Bit::Vector;

## This method gets two files (.newick tree and .fasta with sequences of all its nodes)
## and returns a map of mutations (node name -> array of substitution structs).

sub mutmap_from_files {
	my $tree_file = shift; 
	my $nodeseqs_file = shift;
	
	# read .fasta into a hash
	my %nodeseqs;
	my $seqio = Bio::SeqIO->new(-file => $nodeseqs_file, -format => "fasta");
	while ( my $seqobj = $seqio->next_seq ) {
		my $trimmed_id = (split(/\//, $seqobj->display_id))[0];
    	$nodeseqs{ $trimmed_id } = $seqobj->seq;
	}
	
	# read .newick into Bio::Phylo::IO
	my $tree = parse_tree ($tree_file);
	
	return mutmap($tree, \%nodeseqs);
}


# read .fasta into a hash
sub parse_fasta {
	my $nodeseqs_file = shift;
	my %nodeseqs;
	my $seqio = Bio::SeqIO->new(-file => $nodeseqs_file, -format => "fasta");
	while ( my $seqobj = $seqio->next_seq ) {
		my $trimmed_id = (split(/\//, $seqobj->display_id))[0];
    	$nodeseqs{ $trimmed_id } = $seqobj->seq;
	}
}

sub mutmap {

	my $tree = $_[0];
	my %nodeseqs = %{$_[1]};
	my %map;
	
	my @nodes = $tree -> get_nodes;
	foreach my $node(@nodes){
		if ($node->is_root()) {next;}
		my $name = $node ->get_name();
		#my @nsyn = nsyn_substitutions($nodeseqs{$node->get_ancestors()->[0]->get_name()},
		#							  $nodeseqs{$name});
		my %nsyn = nsyn_substitutions($nodeseqs{$node->get_ancestors()->[0]->get_name()},
									  $nodeseqs{$name});
		$map{$name}=\%nsyn;
	}

	return %map;
};


sub parse_tree {
					my $tree_file = $_[0];
					open TREE, "< $tree_file" or die "Cannot open file ".$tree_file."\n";
					# get a newick string from some source
 					my $tree_string;
 					my $t;
 					while($t = <TREE>){
 						$t =~ s/\n$//;
 						$tree_string = $tree_string.$t;
 					}
 					 close TREE;

 					# Call class method parse from Bio::Phylo::IO
 					# note: newick parser returns 'Bio::Phylo::Forest'
                    # Call ->first to retrieve the first tree of the forest.
 					my $tree = Bio::Phylo::IO->parse(
  					  -string => $tree_string,
   					  -format => 'newick'
 					)->first;

 					return $tree;
	
} 

sub test {
my %mm = mutmap_from_files("C:/Users/weidewind/Documents/CMD/Coevolution/Influenza/Kryazhimsky11/h1.l.r.newick","C:/Users/weidewind/Documents/CMD/Coevolution/Influenza/Kryazhimsky11/h1.all.fa");  
print Dumper ($mm{"INTNODE2340"}[0]);
print ($mm{"INTNODE2340"}[0]->{"Substitution::position"});
}


sub testv {
my $vec1 = Bit::Vector->new(10);
$vec1->Bit_On(3);
my $string = $vec1->to_Bin();
print "'$string'\n";
my $vect = Bit::Vector->new(11);
$vect -> Interval_Copy($vec1,0,0,10);
$vect -> Bit_On(10);
$string = $vect->to_Bin();
print "'$string'\n";
$vec1 = Bit::Vector->new(0);
print ($vec1 ->Size());
$vect =  Bit::Vector->new($vec1 ->Size()+1);
print ($vect ->Size());
}

my $tree = parse_tree("C:/Users/weidewind/Documents/CMD/Coevolution/Influenza/Kryazhimsky11/h1.l.r.newick");
my $fasta = parse_fasta("C:/Users/weidewind/Documents/CMD/Coevolution/Influenza/Kryazhimsky11/h1.all.fa");

my %mutmap = mutmap($tree, $fasta);
compute_bitvectors($tree, \%mutmap, 1);


sub compute_bitvectors{
	my($tree, %mutmap, $site_index) = @_;
	my $root = $tree->get_root();
	$root->set_generic("-mutations_on_path" => Bit::Vector->new(0));
	$root->visit_breadth_first(
				-in => sub{
					my $node=shift;
					my $pnode = $node->get_parent();
					my $pbit_vect = $pnode -> get_generic("-mutations_on_path");
					my $plength = $pbit_vect->Size();
					my $bit_vect = Bit::Vector->new($plength+1);
					$bit_vect -> Interval_Copy($pbit_vect,0,0,$plength);
					if (exists $mutmap{$node->get_name()}{$site_index}) $bit_vect -> Bit_On($plength);
					$node -> set_generic("-mutations_on_path" => $bit_vect);
				},
				-post => sub{ #all daughters have been processed
					my $node=shift;
				}
			);
}


