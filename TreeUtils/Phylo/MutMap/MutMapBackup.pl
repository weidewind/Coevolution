#!/usr/bin/perl
## This class provides a map of mutations (node_name -> array of sites that mutated on the branch),  
## given a tree and sequences of all nodes (including internal ones)

use strict;
use Bio::Phylo::IO;
use DnaUtilities::compare qw(nsyn_substitutions syn_substitutions nsyn_substitutions_codons is_neighbour_changing);
use TreeUtils::Phylo::PhyloUtils qw(remove_zero_branches);

use Bio::Tools::CodonTable;
use TreeUtils::Phylo::FigTree;
use Bio::SeqIO;
use Data::Dumper;
use Bit::Vector;
use Statistics::Test::WilcoxonRankSum;
use Try::Tiny;
use List::Util qw(sum);
use Const::Fast;
	use Switch;

use List::Util qw/shuffle/; 
use Statistics::Basic qw(:all);
use Statistics::TTest;
use Statistics::Descriptive;
use Storable;

use Class::Struct;

$| = 1;

my $static_tree;
my %static_hash_of_nodes;
my %static_fasta;
my $static_protein;
	my %static_subs_on_node;
	my %static_nodes_with_sub;
	
	my %static_background_subs_on_node; # syn for nsyn and vice versa
	my %static_background_nodes_with_sub;
	
	my $static_state; #s=syn or n=nsyn
	
	my %static_ring_hash; 
	my %static_depth_hash; 
	my $static_alignment_length;


	my @static_sorted_nodnames;
	my @static_sorted_sites;
	
	my %static_subtree_info;
	my %static_distance_hash;
	
{
	my %distance_hash;
	sub set_node_distance {
		$distance_hash{$_[0]}->{$_[1]} = $_[2];	
	}
	
	sub set_alignment_length {
		$static_alignment_length = $_[0]; 
	}
	sub has_node_distance {
		if (!defined $distance_hash{$_[0]}->{$_[1]}){
			return 0;
		}
		else {
			return 1;
		}
	}
	sub get_node_distance {
		return $distance_hash{$_[0]}->{$_[1]};
	}
	
	
struct Mutation => {
	site_index => '$',
	node => '$',
};
	
	sub set_mutmap {
		$static_tree = parse_tree("C:/Users/weidewind/Documents/CMD/Coevolution/Influenza/Kryazhimsky11/".$_[0].".l.r.newick");
		my @arr = parse_fasta("C:/Users/weidewind/Documents/CMD/Coevolution/Influenza/Kryazhimsky11/".$_[0].".all.fa");
		my %fasta = %{$arr[0]};
		my $alignment_length = $arr[1];
		$static_alignment_length = $alignment_length;
		$static_protein  = $_[0];
		%static_fasta = %fasta;
			
		my @nodes = $static_tree -> get_nodes;
		foreach my $node(@nodes){
			#if ($node->is_root()) {next;}
			my $name = $node ->get_name();
			$static_hash_of_nodes{$name} = \$node;
		}
		
		my @mutmaps;
		my @bkg_mutmaps;
		if($_[1] eq "syn"){
			 @mutmaps = synmutmap($static_tree, \%fasta);
			 @bkg_mutmaps = codonmutmap($static_tree, \%fasta);
			 $static_state = "s";
		} 
		elsif($_[1] eq "nsyn"){
			 @mutmaps = codonmutmap($static_tree, \%fasta);
			 @bkg_mutmaps = synmutmap($static_tree, \%fasta);
			 $static_state = "n";
		} 
		else {
			die "only syn or nsyn can be used as the second argument; unknown $_[1] was used instead";
		}
		%static_subs_on_node = %{$mutmaps[0]};
		%static_nodes_with_sub = %{$mutmaps[1]};
		
		%static_background_subs_on_node = %{$bkg_mutmaps[0]};
		%static_background_nodes_with_sub = %{$bkg_mutmaps[1]};
	
		
		
	}
	
	sub get_static_tree{
		return  $static_tree;
	}
	
	sub get_static_subs_on_node{
		return  %static_subs_on_node;
	}
	
	sub get_static_nodes_with_sub{
		return  %static_nodes_with_sub;
	}
	

}



## according to http://www.biomedcentral.com/1471-2148/10/253
		const my @n1_decreasing => ("AGG", "TCG", "GAT", "CGT", "ACC", "GCC", "CAG", "GGG", "GGC");

		const my @h1_decreasing => ("ACG", "TCA", "CTC", "GCG", "GCA", "CCG", "TGC", "GTG");

		const my @n2_decreasing => ("AAT", "CTC", "GAG", "TCT", "ACT", "TGT", "CCG", "GGC", "GAC", "AAA", "TCA");

		const my @h3_decreasing => ("CTG", "CGC", "CCT", "TGC",  "GAC", "AGG", "TAT", "AAG", "GGG", "CGG");

		const my @n1_increasing => ("AGA", "ACA", "GGA", "CAC", "TCA", "CTT", "CAA", "AGT");

		const my @h1_increasing  => ("ACA", "GCC", "CCT", "TGT", "TCC", "AGC");

		const my @n2_increasing => ("AAC", "TCC", "GAA", "GTT", "TGC", "GAT", "AAG", "GCC", "ACA");
		
		const my @h3_increasing => ("TTG", "AGA", "TGT", "GCC",  "CTA", "GAT", "TAC", "CCG", "GGA", "AAA", "CCC");

		const my @all_codons => ("TCA", "TCC", "TCG", "TCT", "TTC", "TTT", "TTA", "TTG", "TAC", "TAT", "TAA", "TAG", "TGC", "TGT", "TGA", "TGG", "CTA", "CTC", "CTG", "CTT", "CCA", "CAT", "CAA", "CAG", "CGA", "CGC", "CGG", "CGT", "ATA", "ATC", "ATT", "ATG", "ACA", "ACC", "ACG", "ACT", "AAC", "AAT", "AAA", "AAG", "AGC", "AGT", "AGA", "AGG", "CCC", "CCG", "CCT", "CAC", "GTA", "GTC", "GTG", "GTT", "GCA", "GCC", "GCG", "GCT", "GAC", "GAT", "GAA", "GAG", "GGA", "GGC", "GGG", "GGT");

		

## returns a hash: key - codon, value - -1, if it is decreasing over time,  
##                                       1, if it is increasing, 
##										 0 otherwise.
sub codon_evolution{
	my $protein = $_[0];
	my %hash;
	@hash{@all_codons} = 0;
	switch($protein){
		case "n1" {	@hash{@n1_decreasing} = -1; 
			        @hash{@n1_increasing} = 1; }
		case "n2" {	@hash{@n2_decreasing} = -1; 
					@hash{@n2_increasing} = 1; }
		case "h1" {	@hash{@h1_decreasing} = -1; 
					@hash{@h1_increasing} = 1; }
		case "h3" {	@hash{@h3_decreasing} = -1; 
					@hash{@h3_increasing} = 1; }
		}
	
	return %hash;
}



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
	my $length;
	while ( my $seqobj = $seqio->next_seq ) {
		my $trimmed_id = (split(/\//, $seqobj->display_id))[0];
    	$nodeseqs{ $trimmed_id } = $seqobj->seq;
    	if (!$length){
    		$length = $seqobj->length();
    	}
	}
	return (\%nodeseqs, $length);
}

sub mutmap {

	my $tree = $_[0];
	my %nodeseqs = %{$_[1]};
	my %subs_on_node;
	my %nodes_with_sub;
	
	my @nodes = $tree -> get_nodes;
	foreach my $node(@nodes){
		if ($node->is_root()) {next;}
		my $name = $node ->get_name();
		#my @nsyn = nsyn_substitutions($nodeseqs{$node->get_ancestors()->[0]->get_name()},
		#							  $nodeseqs{$name});
		my %nsyn = nsyn_substitutions($nodeseqs{$node->get_ancestors()->[0]->get_name()},
									  $nodeseqs{$name});					  
		$subs_on_node{$name}=\%nsyn;
		for	my $site_index(keys %nsyn){
			if (! exists $nodes_with_sub{$site_index}){
				$nodes_with_sub{$site_index} = ();
			}
			push (@{$nodes_with_sub{$site_index}}, \$node);
		}	
	}

	return (\%subs_on_node, \%nodes_with_sub);
};

sub codonmutmap {

	my $tree = $_[0];
	my %nodeseqs = %{$_[1]};
	my %subs_on_node;
	my %nodes_with_sub;
	
	my @nodes = $tree -> get_nodes;
	foreach my $node(@nodes){
		if ($node->is_root()) {next;}
		my $name = $node ->get_name();
		#my @nsyn = nsyn_substitutions($nodeseqs{$node->get_ancestors()->[0]->get_name()},
		#							  $nodeseqs{$name});
		my %nsyn = nsyn_substitutions_codons($nodeseqs{$node->get_ancestors()->[0]->get_name()},
									  $nodeseqs{$name});					  
		$subs_on_node{$name}=\%nsyn;
		for	my $site_index(keys %nsyn){
			if (! exists $nodes_with_sub{$site_index}){
				$nodes_with_sub{$site_index} = ();
			}
			push (@{$nodes_with_sub{$site_index}}, \$node);
		}	
	}

	return (\%subs_on_node, \%nodes_with_sub);
};


sub synmutmap {

	my $tree = $_[0];
	my %nodeseqs = %{$_[1]};
	my %subs_on_node;
	my %nodes_with_sub;
	
	my @nodes = $tree -> get_nodes;
	foreach my $node(@nodes){
		if ($node->is_root()) {next;}
		my $name = $node ->get_name();
		my %nsyn = syn_substitutions($nodeseqs{$node->get_ancestors()->[0]->get_name()},
									  $nodeseqs{$name});					  
		$subs_on_node{$name}=\%nsyn;
		for	my $site_index(keys %nsyn){
			if (! exists $nodes_with_sub{$site_index}){
				$nodes_with_sub{$site_index} = ();
			}
			push (@{$nodes_with_sub{$site_index}}, \$node);
		}	
	}

	return (\%subs_on_node, \%nodes_with_sub);
}

sub incidence_matrix_1 {

	my %subs_on_node = %static_subs_on_node;
	my %matrix;
	my $length;
	if ($static_state eq 'n'){
		$length = $static_alignment_length/3;
	}
	else {
		$length = $static_alignment_length;
	}

	my @nodes = $static_tree -> get_nodes;
	foreach my $node(@nodes){
		my $name = $node ->get_name();
		my $incidence_vector =  Bit::Vector->new($length+1);
		
		foreach my $index(keys %{$subs_on_node{$name}}){
			print "index ".$index;
			$incidence_vector-> Bit_On($index); #indices start from 1
		}
		$matrix{$name} = $incidence_vector;
		print $matrix{$name}->to_Bin();
		print "\n";
	}
	
};

sub complete_incidence_matrix {
	my %subs_on_node = %static_subs_on_node;
	my %nodes_with_sub = %static_nodes_with_sub;
	my %matrix;
	my $length;
	if ($static_state eq "n"){
		$length = $static_alignment_length/3;
	}
	elsif ($static_state eq "s") {
		$length = $static_alignment_length;
	}
	
	my @sorted_sites;
	my @sorted_nodnames;
	
	# select nodes with at least one mutation of the corresponding type (syn or nsyn, depending on the mutmap state)
	my @nodes = $static_tree -> get_nodes;
	foreach my $node(@nodes){
		my $name = $node ->get_name();
		
			push @sorted_nodnames, $name;
		
	}
	
	my %empty_nodes_hash = map { $_ => 0 } @sorted_nodnames;
	my %incidence_hash;
	
	# select sites with at least 3 mutations of the corresponding type
	# upd - not sure if such sites should be excluded, changed minimum to 1
	foreach my $ind(1..$length){
			my %site_incidence = %empty_nodes_hash;
			push @sorted_sites, $ind;
			#print " added $ind to sorted sites\n";
			foreach my $node(@{$nodes_with_sub{$ind}}){
				$site_incidence{$$node->get_name()} = 1;
			}
			$incidence_hash{$ind} = \%site_incidence;
		
	}
	
	@static_sorted_nodnames = @sorted_nodnames;
	@static_sorted_sites = @sorted_sites;
	
	return %incidence_hash;
}


sub incidence_matrix {

	my %subs_on_node = %static_subs_on_node;
	my %nodes_with_sub = %static_nodes_with_sub;
	my %matrix;
	my $length;
	if ($static_state eq "n"){
		$length = $static_alignment_length/3;
	}
	elsif ($static_state eq "s") {
		$length = $static_alignment_length;
	}
	
	my @sorted_sites;
	my @sorted_nodnames;
	
	# select nodes with at least one mutation of the corresponding type (syn or nsyn, depending on the mutmap state)
	my @nodes = $static_tree -> get_nodes;
	foreach my $node(@nodes){
		my $name = $node ->get_name();
		if (scalar (keys %{$subs_on_node{$name}}) > 0){
			push @sorted_nodnames, $name;
		}
	}
	
	my %empty_nodes_hash = map { $_ => 0 } @sorted_nodnames;
	my %incidence_hash;
	
	# select sites with at least 3 mutations of the corresponding type
	# upd - not sure if such sites should be excluded, changed minimum to 1
	foreach my $ind(1..$length){
		if($nodes_with_sub{$ind} && scalar @{$nodes_with_sub{$ind}} > 0){# possibility to have at least 2 mutations after an ancestral one
			my %site_incidence = %empty_nodes_hash;
			push @sorted_sites, $ind;
			#print " added $ind to sorted sites\n";
			foreach my $node(@{$nodes_with_sub{$ind}}){
				$site_incidence{$$node->get_name()} = 1;
			}
			$incidence_hash{$ind} = \%site_incidence;
		}
	}
	
	@static_sorted_nodnames = @sorted_nodnames;
	@static_sorted_sites = @sorted_sites;
	
	return %incidence_hash;
};

sub print_incidence_matrix {
	my %incidence_hash = %{$_[0]};
	my $path = $_[1];
	my $matrix_file = $path.$static_protein."_incidence_matrix";
	my $sorted_sites_file = $path.$static_protein."_sorted_sites";
	my $sorted_nodnames_file = $path.$static_protein."_sorted_nodnames";
	open MATRIX, ">$matrix_file" or die "Cannot open file ".$matrix_file."\n";
	foreach my $nodname (@static_sorted_nodnames){
		foreach my $ind (@static_sorted_sites){
			print MATRIX $incidence_hash{$ind}->{$nodname};
		}
		print MATRIX "\n";
	}
	close MATRIX;
	
	open SSITES, ">$sorted_sites_file" or die "Cannot open file ".$sorted_sites_file."\n";
	foreach my $ind(@static_sorted_sites){
		print SSITES $ind."\n";
	}
	close SSITES;
	
	open SNODES, ">$sorted_nodnames_file" or die "Cannot open file ".$sorted_nodnames_file."\n";
	foreach my $name(@static_sorted_nodnames){
		print SNODES $name."\n";
	}
	close SNODES;
	
}

sub read_incidence_matrix {
	my $matrix_file = $_[0];
	open MATRIX, "<$matrix_file" or die "Cannot open file ".$matrix_file."\n";
	my %subs_on_node;
	my %nodes_with_sub;
	my $line_index = 0;
	while(<MATRIX>){
		if (/^$/) {last;}
			my $nodname = $static_sorted_nodnames[$line_index];
			my @sites = split(',');
			my %substs;
			foreach my $s(@sites){
				my $ind = $static_sorted_sites[$s-1];
#				print " $s is $ind\n";
				my $p=Substitution->new();
				$p->position($ind);
				$p->ancestral_allele("ATG");
				$p->derived_allele("ATG");
				$substs{$ind} = $p;
				if (! exists $nodes_with_sub{$ind}){
					$nodes_with_sub{$ind} = ();
				}
	#			print "\n nodndame ".$_."\n";
	#			print "\nREF 1 ".ref($static_hash_of_nodes{$nodname})."\n";
	#			print "\nREF 2 ".ref(${$static_hash_of_nodes{$nodname}})."\n";
				push (@{$nodes_with_sub{$ind}}, \${$static_hash_of_nodes{$nodname}}); #�������� �� ������ �� �����
			}
			$subs_on_node{$nodname} = \%substs;
			$line_index++;

	}
	close MATRIX;
	return (\%subs_on_node, \%nodes_with_sub);
			
}

sub test_mutmaps{
	
	foreach my $ind(keys %static_nodes_with_sub){
		my $nodes_count = scalar @{$static_nodes_with_sub{$ind}};
		print " index $ind, $nodes_count nodes ";
		foreach my $node(@{$static_nodes_with_sub{$ind}}){
			print $$node->get_name()."\t";
		}
		print "\n";
	}
	
	foreach my $nodname(keys %static_subs_on_node){
		my $sites_count = scalar keys $static_subs_on_node{$nodname};
		print " node $nodname, $sites_count sites\n";
	}
}

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

sub test2 {
my %mm = mutmap_from_files("C:/Users/weidewind/Documents/CMD/Coevolution/Influenza/Kryazhimsky11/h1.l.r.newick","C:/Users/weidewind/Documents/CMD/Coevolution/Influenza/Kryazhimsky11/h1.all.fa");  
print Dumper ($mm{"INTNODE2340"});
#print ($mm{"INTNODE2340"}[0]->{"Substitution::position"});
}


sub testv {
my $vec1 = Bit::Vector->new(10);
$vec1->Bit_On(3);
my $string = $vec1->to_Bin();
print "'$string'\n";
my $vect = Bit::Vector->new(11);
$vect -> Interval_Copy($vec1,0,3,10);
$vect -> Bit_On(10);
$string = $vect->to_Bin();
print "'$string'\n";
$vect->Move_Left(1);
$string = $vect->to_Bin();
print "'$string'\n";
my $vect2 = Bit::Vector->new(6);
$vect2 -> Bit_On(3);
my $vect3 = Bit::Vector->new(6);
$vect3 -> Bit_On(2);
my $vect4 = Bit::Vector->new(6);
$vect4->Divide($vect2, $vect3, Bit::Vector->new(6));
print ($vect2->to_Bin()."\t".$vect3->to_Bin()."\t".$vect4->to_Bin()."\n");

}

sub wilkotest{
	my @a1;
	my $undefi;
	push @a1, (1,2,7,9);
	#push @a1, $undefi;
	push @a1, 3;
	my @a2;
	push @a2, (5,2,8,13);
	my $wilcox_test = Statistics::Test::WilcoxonRankSum->new();
	$wilcox_test->load_data(\@a1, \@a2);
  my $prob = $wilcox_test->probability();
  my $pf = sprintf '%f', $prob; # prints 0.091022
  print "\n";
  print $wilcox_test->probability_status();
  print $wilcox_test->summary();
    print "\n";
    #ok, undef values don't influence the statistics
}

sub patrtest{
		my $tree = parse_tree("C:/Users/weidewind/Documents/CMD/Coevolution/Influenza/Kryazhimsky11/n1.l.r.newick");
	foreach my $node(@{$tree->get_internals()}){
		if ($node->get_name() eq "INTNODE2473"){
		foreach my $n(@{$tree->get_internals()}){
			if ($n->get_name() eq "INTNODE2479"){
				print get_mrcn($node, $n)->get_name()."\n";
			print $node->get_name()."\t".$n->get_name()."\t".calc_true_patristic_distance($node, $n)."\n";
			}
	}
		}
	}
}


## unlike original bio::phylo get_mrca, returns the node n1 closest to the root, if you give it two sequential nodes n1 and n2
## (in that case get_mrca returns the youngest ancestor of n1)

sub get_mrcn {
        my ( $node, $other_node ) = @_;
        if ( $node->get_id == $other_node->get_id ) {
            return $node;
        }
        my $self_anc  = $node->get_ancestors;
		unshift @{$self_anc}, $node;
        my $other_anc = $other_node->get_ancestors;
		unshift @{$other_anc}, $other_node;
        for my $i ( 0 .. $#{$self_anc} ) {
            my $self_anc_id = $self_anc->[$i]->get_id;
            for my $j ( 0 .. $#{$other_anc} ) {
                if ( $self_anc_id == $other_anc->[$j]->get_id ) {
                    return $self_anc->[$i];
                }
            }
        }

        return $self_anc->[-1];
    }

## the only difference from calc_patristic_distance is that it uses get_mrcn instead of get_mrca
## If you give it two sequential nodes, it returns the distance between them 

 sub calc_true_patristic_distance {
        my ( $node, $other_node ) = @_;
        my $patristic_distance = 0;
        my $mrca    = get_mrcn($node, $other_node);
        my $mrca_id = $mrca->get_id;
        while ( $node->get_id != $mrca_id ) {
            my $branch_length = $node->get_branch_length;
            if ( defined $branch_length ) {
                $patristic_distance += $branch_length;
            }
            $node = $node->get_parent;
        }
        while ( $other_node and $other_node->get_id != $mrca_id ) {
            my $branch_length = $other_node->get_branch_length;
            if ( defined $branch_length ) {
                $patristic_distance += $branch_length;
            }
            $other_node = $other_node->get_parent;
        }
        return $patristic_distance;
    }
    
    
## like calc_patristic_distance, but returns 0 for two sequential nodes    
    sub calc_my_distance {
        my ( $node, $other_node ) = @_;
        my $patristic_distance = 0;
        my $mrca    = get_mrcn($node, $other_node);
        my $mrca_id = $mrca->get_id;
        if ( $node->get_id == $mrca_id || $other_node->get_id == $mrca_id){
        	return 0;
        }
        while ( $node->get_id != $mrca_id ) {
            my $branch_length = $node->get_branch_length;
            if ( defined $branch_length ) {
                $patristic_distance += $branch_length;
            }
            $node = $node->get_parent;
        }
        while ( $other_node and $other_node->get_id != $mrca_id ) {
            my $branch_length = $other_node->get_branch_length;
            if ( defined $branch_length ) {
                $patristic_distance += $branch_length;
            }
            $other_node = $other_node->get_parent;
        }
        return $patristic_distance;
    }
    
    sub node_distance {
    	my ( $node, $other_node ) = @_;
    	if  (has_node_distance($node, $other_node)){
    		return get_node_distance($node, $other_node);
    	}
    	else {
    		## calc_true instead of calc_my since 02 06 2015
    		my $dist = calc_true_patristic_distance($node, $other_node);
    		set_node_distance($node, $other_node, $dist);
    		return $dist;
    	}
    	
    }


#set_mutmap("h3", "nsyn");
#my @sites = (229,291,469);
#for my $site(@sites){
#	print_static_tree_with_mutations($site);
#}
#print_tree_with_mutations(503);
#print_tree_with_mutations(465);
# prints nexus tree, on wich all mutations in the specified site are shown 

sub print_static_tree_with_mutations{
my $site = shift;
my $myCodonTable   = Bio::Tools::CodonTable->new();

my %sites;
my %color;
foreach my $n(@{$static_nodes_with_sub{$site}}){
	my $sub = ${$static_subs_on_node{$$n->get_name()}}{$site};
	$sites{$$n->get_name()} = $$n->get_name()."_".$site."_".$sub->{"Substitution::ancestral_allele"}.
							  "(".$myCodonTable->translate($sub->{"Substitution::ancestral_allele"}).")->".
							  $sub->{"Substitution::derived_allele"}.
							  "(".$myCodonTable->translate($sub->{"Substitution::derived_allele"}).")";
	$color{$$n->get_name()} = "-16776961";
}

open TREE, ">".$static_protein."_sites_".$site.".tre";
print TREE "#NEXUS\n\nbegin trees;\n";
print TREE "\ttree $site = [&R] ";
my $tree_name=tree2str($static_tree,sites => \%sites, color=>\%color);
print TREE $tree_name;
print TREE "\nend;";
close TREE;

foreach my $trn(@{$static_nodes_with_sub{$site}}){
	print ($$trn->get_name()."\t");
	print (${$static_subs_on_node{$$trn->get_name()}}{$site}->{"Substitution::derived_allele"});
	print "\n";
	foreach my $trr(@{$static_nodes_with_sub{$site}}){
		
		print "\t".calc_true_patristic_distance($$trr, $$trn)."_";
		print (${$static_subs_on_node{$$trr->get_name()}}{$site}->{"Substitution::derived_allele"}."_");
		print $$trr->get_name();
		print "\n";
	}
	print "\n";
}
}


sub print_tree_with_mutations{
my $site = shift;
my $prot = shift;
my $tree = parse_tree("C:/Users/weidewind/Documents/CMD/Coevolution/Influenza/Kryazhimsky11/$prot.l.r.newick");
my %fasta = parse_fasta("C:/Users/weidewind/Documents/CMD/Coevolution/Influenza/Kryazhimsky11/$prot.all.fa");
my @mutmaps = codonmutmap($tree, \%fasta);
my $myCodonTable   = Bio::Tools::CodonTable->new();
my %subs_on_node = %{$mutmaps[0]};
my %nodes_with_sub = %{$mutmaps[1]};

foreach my $s( keys %nodes_with_sub){
	print $s."\t";
}

my %sites;
my %color;
foreach my $n(@{$nodes_with_sub{$site}}){
	my $sub = ${$subs_on_node{$$n->get_name()}}{$site};
	$sites{$$n->get_name()} = $$n->get_name()."_".$site."_".$sub->{"Substitution::ancestral_allele"}.
							  "(".$myCodonTable->translate($sub->{"Substitution::ancestral_allele"}).")->".
							  $sub->{"Substitution::derived_allele"}.
							  "(".$myCodonTable->translate($sub->{"Substitution::derived_allele"}).")";
	$color{$$n->get_name()} = "-16776961";
}

open TREE, ">".$prot."_sites_".$site.".tre";
print TREE "#NEXUS\n\nbegin trees;\n";
print TREE "\ttree $site = [&R] ";
my $tree_name=tree2str($tree,sites => \%sites, color=>\%color);
print TREE $tree_name;
print TREE "\nend;";
close TREE;

foreach my $trn(@{$nodes_with_sub{$site}}){
	print ($$trn->get_name()."\t");
	print (${$subs_on_node{$$trn->get_name()}}{$site}->{"Substitution::derived_allele"});
	print "\n";
	foreach my $trr(@{$nodes_with_sub{$site}}){
		
		print "\t".calc_true_patristic_distance($$trr, $$trn)."_";
		print (${$subs_on_node{$$trr->get_name()}}{$site}->{"Substitution::derived_allele"}."_");
		print $$trr->get_name();
		print "\n";
	}
	print "\n";
}
}

sub mean_ignore_nulls{
	if (!defined $_[0]){
		return 0;
	}

	my @arr = @{$_[0]};
	my $count = 0;
	my $sum = 0;
	for my $num(@arr){
		if (defined $num){
			$count++;
			$sum += $num;
		}
	}
	if ($count == 0){
		return 0;
	}
	return $sum/$count;
	
}


my @h3_host_shift = (2,	3,	4,	9,	10,	11,	14,	16,	18,	19,	20,	22,	23,	25,	47,	66,	69,	73,	78,	79,	83,	97,	98,	99,	108,	110,	137,	142,	153,	159,	160,	161,	162,	176,	179,	206,	208,	212,	229,	230,	238,	244,	260,	264,	285,	291,	323,	328,	329,	400,	402,	462,	492,	506,	541);
my @h1_host_shift = (2, 9, 14, 15, 22, 47, 61, 62, 71, 73, 78, 85, 88, 89, 97, 100, 102, 113, 130, 132, 138, 144, 146, 149, 151, 153, 154, 155, 157, 167, 168, 169, 171, 172, 173, 176, 177, 180, 186, 200, 201, 202, 203, 205, 206, 209, 210, 211, 212, 218, 223, 227, 230, 232, 235, 238, 240, 243, 251, 252, 257, 261, 268, 274, 275, 277, 278, 285, 286, 288, 289, 293, 294, 299, 302, 314, 323, 324, 325, 326, 331, 337, 389, 415, 420, 434, 453, 459, 466, 470, 515, 541, 542, 548);
my @n1_host_shift = (3, 5, 8, 12, 13, 14, 16, 20, 26, 29, 34, 40, 41, 42, 43, 46, 47, 51, 52, 53, 59, 64, 66, 67, 69, 70, 71, 72, 74, 75, 76, 78, 79, 80, 81, 82, 83, 85, 93, 95, 99, 101, 105, 111, 114, 116, 136, 149, 157, 189, 195, 200, 206, 210, 211, 214, 220, 221, 222, 223, 232, 241, 250, 257, 258, 263, 264, 267, 273, 274, 285, 287, 288, 289, 309, 311, 329, 339, 340, 341, 351, 354, 355, 365, 367, 369, 382, 386, 388, 390, 393, 394, 396, 427, 430, 432, 434, 451, 454, 455);
my @n2_host_shift = (7, 9, 19, 22, 24, 26, 28, 31, 33, 38, 39, 40, 41, 42, 44, 45, 48, 50, 51, 52, 57, 58, 59, 60, 62, 66, 69, 70, 72, 73, 77, 79, 81, 83, 85, 86, 93, 95, 100, 113, 116, 125, 126, 143, 147, 149, 150, 155, 187, 192, 199, 206, 210, 212, 216, 220, 221, 234, 238, 257, 267, 275, 283, 284, 286, 290, 296, 305, 308, 310, 311, 312, 313, 315, 328, 331, 332, 336, 338, 342, 347, 356, 360, 367, 368, 369, 370, 378, 380, 381, 384, 385, 386, 390, 393, 396, 399, 400, 401, 403, 415, 431, 435, 437, 445, 466);
my @h1_host_shift_001 = (203, 168, 299, 251, 288, 201, 167, 252, 302, 62, 9, 238, 314, 324, 275, 285, 154, 172, 176, 459, 420, 2, 211, 202, 130, 470, 274, 257, 14, 323, 89, 294, 261, 235, 100, 286, 415, 200, 206, 15, 85, 78, 210, 71, 453, 466, 337, 22);
my @h3_host_shift_001 = (16, 108, 229, 244, 79, 73, 83, 161, 260, 20, 9 );
my @n2_host_shift_001 = (386, 384, 381, 328, 83, 70, 81, 192, 51, 147, 125, 283, 41, 286, 77, 72, 378, 331, 126, 155, 50, 62, 338, 369, 60, 315, 216, 399, 396) ;
my @n1_host_shift_001 = (189, 382, 214, 340, 311, 274, 157, 430, 455, 74, 341, 220, 221, 288, 351, 80, 264, 289, 365, 339, 52, 46, 59, 42, 34, 47, 393, 427, 3, 67, 309, 329, 29);

my @h1_epitopes = qw(141 142 171 173 175 176 178 179 180 169 172 205 206 209 211 153 156 158 182 186 195 220 237 238 253 286 87 88 90 91 92 132);
my @h1_increased_binding = qw(141 142 171 178 179 180 169 193 205 209 211 156 237 87 88 132 257 175 158 106);
my @n1_epitopes = qw(380 381 382 383 384 385 386 388 389 390 393 397 398 199 200 201 202 223 329 330 332 331 333 336 337 339 340 341 343 344 356 363 364 365 366 367);
my @h3_epitopes = qw( 138 140 142 147 148 146 149 151 153 154 156 158 159 160 161 162 166 168 184 144 145 171 172 173 174 175 176 179 180 181 202 203 204 205 206 208 209 210 212 213 214 60 61 62 63 64 66 67 69 70 289 291 292 294 295 296 310 313 315 316 320 321 323 324 325 326 327 328 112 118 119 133 137 183 186 187 188 189 190 191 192 193 195 198 217 219 223 224 225 228 229 230 231 232 233 234 235 242 243 244 245 246 254 256 258 260 262 263 264 73 75 78 79 83 91 94 96 97 98 99 102 103 104 107 108 110 125 276 277 278 281 );
my @n2_epitopes = qw(383 384 385 386 387 389 390 391 392 393 394 396 399 400 401 403 197 198 199 200 221 222 328 329 330 331 332 334 336 338 339 341 342 343 344 346 347 357 358 359  366 367 368 369 370);
my @n1_wan_epitopes = qw(248, 249, 250, 273, 309, 338, 339, 341, 343, 396, 397, 456);
#h1 Huang (antigenic), as is in file Tables_main (Huang + 17, from msa)
my @h1_antigenic = qw( 138 144 145 147 150 158 163 142 170 177 200 203 206 207 208 210 211 52 53 60 288 290 291 294 312 327 111 180 222 226 233 239 241 64 71 86 88 90 97 99 284 );


#h1 Ren, Li, Liu (antigenic) - intersection of two methods
my @h1_antigenic_ren = qw(60 71	88	138	142	144	147	158	204	207	210	222	338);

my @h1_antigenic_Huang_and_host_shift = qw(71	200	203	206	210	211	288	294	);




my @n2_surface = (34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 88, 89, 90, 91, 92, 93, 95, 107, 110, 111, 112, 113, 118, 125, 126, 127, 128, 130, 141, 143, 146, 147, 149, 150, 151, 152, 153, 154, 160, 161, 162, 169, 171, 173, 187, 189, 196, 197, 198, 199, 200, 208, 209, 210, 212, 215, 216, 218, 219, 220, 221, 222, 224, 234, 236, 244, 245, 246, 247, 248, 249, 250, 251, 253, 258, 259, 261, 262, 263, 264, 265, 267, 268, 269, 270, 271, 273, 277, 283, 284, 285, 286, 292, 295, 296, 304, 306, 307, 308, 309, 310, 311, 312, 313, 315, 326, 328, 329, 330, 331, 332, 334, 336, 337, 338, 339, 341, 342, 343, 344, 346, 347, 356, 357, 358, 359, 366, 367, 368, 369, 370, 371, 378, 380, 381, 383, 384, 385, 386, 387, 388, 389, 390, 391, 392, 393, 394, 396, 399, 400, 401, 402, 403, 413, 414, 415, 416, 417, 430, 431, 432, 433, 434, 435, 437, 450, 451, 452, 453, 455, 456, 457, 459, 461, 463, 464, 465, 466, 468, 469, 470);
my @n1_surface = (34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 86, 88, 89, 90, 93, 94, 95, 111, 118, 126, 127, 128, 136, 141, 143, 146, 147, 148, 149, 150, 151, 152, 154, 162, 163, 165, 172, 174, 189, 191, 199, 200, 201, 202, 209, 210, 211, 214, 215, 217, 218, 220, 221, 222, 223, 224, 226, 236, 237, 247, 248, 249, 250, 251, 252, 255, 258, 260, 261, 263, 264, 265, 266, 267, 268, 269, 271, 273, 274, 275, 279, 285, 286, 287, 288, 290, 296, 297, 298, 304, 306, 307, 308, 309, 311, 312, 313, 314, 326, 328, 329, 330, 331, 332, 333, 334, 335, , , 336, 337, 338, 340, 341, 349, 350, 351, 352, 360, 361, 362, 363, 364, 365, 372, 374, 375, 378, 379, 380, 381, 382, 383, 384, 385, 386, 387, 389, 391, 392, 393, 394, 395, 398, 406, 407, 408, 413, 414, 415, 416, 417, 426, 427, 429, 430, 431, 432, 433, 434, 435, 437, 439, 450, 451, 452, 454, 455, 456, 457, 461, 463, 464, 465, 467, 468);

my @n1_internal = (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,85,87,91,92,96,97,98,99,100,101,102,103,104,105,106,107,108,109,110,112,113,114,115,116,117,119,120,121,122,123,124,125,129,130,131,132,133,134,135,137,138,139,140,142,144,145,153,155,156,157,158,159,160,161,164,166,167,168,169,170,172,174,175,176,177,178,179,180,181,182,183,184,185,186,187,189,191,192,193,194,195,196,197,202,203,204,205,206,207,211,212,215,218,224,226,227,228,229,230,231,232,233,234,237,238,239,240,241,242,243,244,245,252,253,255,256,258,261,269,271,275,276,277,279,280,281,282,283,288,290,291,292,293,294,298,299,300,301,302,304,306,310,315,316,317,318,319,320,321,322,323,324,325,327,333,342,345,346,347,348,349,350,351,356,357,358,359,360,361,362,369,370,371,372,373,374,376,379,380,391,393,394,400,401,403,404,405,406,407,408,409,418,419,420,421,422,423,424,425,428,434,439,441,442,443,444,445,446,447,448,449,450,454,459,460,461,463,467,470,471,472);
my @n2_internal = ( 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 87, 94, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 108, 109, 114, 115, 116, 117, 119, 120, 121, 122, 123, 124, 129, 131, 132, 133, 134, 135, 136, 137, 138, 139, 140, 142, 144, 145, 148, 155, 156, 157, 158, 159, 163, 164, 165, 166, 167, 168, 170, 172, 174, 175, 176, 177, 178, 179, 180, 181, 182, 183, 184, 185, 186, 188, 190, 191, 192, 193, 194, 195, 201, 202, 203, 204, 205, 206, 207, 211, 213, 214, 217, 223, 225, 226, 227, 228, 229, 230, 231, 232, 233, 235, 237, 238, 239, 240, 241, 242, 243, 252, 254, 255, 256, 257, 260, 266, 272, 274, 275, 276, 278, 279, 280, 281, 282, 287, 288, 289, 290, 291, 293, 294, 297, 298, 299, 300, 301, 302, 303, 305, 314, 316, 317, 318, 319, 320, 321, 322, 323, 324, 325, 327, 333, 335, 340, 345, 348, 349, 350, 351, 352, 353, 354, 355, 360, 361, 362, 363, 364, 365, 372, 373, 374, 375, 376, 377, 379, 382, 395, 397, 398, 404, 405, 406, 407, 408, 409, 410, 411, 412, 418, 419, 420, 421, 422, 423, 424, 425, 426, 427, 428, 429, 436, 438, 439, 440, 441, 442, 443, 444, 445, 446, 447, 448, 449, 454, 458, 460, 462, 467);

#h3 antigentic Steinbruck
my @h3_antigenic  = qw( 138 160 171 223 161 205 233 294 66 153 174 276 140 151 230 278 78 172 212 292 41 91 99 147 202 218 238 241 19 204 69 180 190 209 217 229 246 98 149 159 162 176 213 18 70 188 260 206 242 );

# h3 antigenic change Koel
my @h3_antigenic_koel = (161, 171, 172, 174, 175, 205, 209);
my @h3_antigenic_smith = (138, 160, 153, 161, 149, 159, 162, 140, 147, 171, 204, 180, 205, 209, 174, 172, 176, 213, 175, 206, 212, 294, 66, 70, 223, 190, 217, 229, 233, 246, 188, 260, 292, 98, 276, 78, 91, 99, 69);

#my @h1_pocket_closest = (207,239,166,211,148,111,208,206,234,149,235,196,203,210,151,150,241,204,237,236,240,238,205,209,242,110,197,195,147,202,212,165,152,112,167,233,158,243,244,198,265,159,199,263,201,153,169,262,160);
my @h1_pocket_closest = (207,239,166,111,149,196,203,241,147,168,240,242,110,169,146,197,195,148,206,238,150,204,208,202,167,165,112,243,244,235,198,265,158,199,263,144,205,151,236,159,201,237,152,145,200,143,113,209,245);
my @h3_pocket_closest = (152,241,242,206,238,210,151,114,153,169,165,205,209,42,29,28,211,239,150,207,237,168,113,154,170,115,240,243,212,208,166,43,162,164,161,171,268,175,27,203,204,163,149,155,244,167,172,41,30,148,236,44,156,269);
my @n2_pocket_closest = (151,277,406,150,407,276,278,152,405,226,350,425,291,225,292);
my @n2_bigger_pocket_closest = (292,152,119,276,224,227,406,118,277,371,117,225,407,275,278,372,370,291,223,153,151,228,405,293,226,120,350,425,180,404,181,242,178,133,134,300,241,440,179,441,365,222,150,349);
my @n1_pocket_closest = (402,151,278,152,403,401,277,279,150,227,425,347,292,424,226,228);

my @h1_surface = (28,30,38,39,40,42,52,53,54,56,62,63,64,68,71,83,86,91,92,100,101,103,111,132,136,137,138,141,142,146,148,150,155,156,157,158,171,172,173,176,178,182,184,186,200,201,202,205,206,211,212,221,227,232,235,237,238,252,253,277,278,283,285,287,289,290,292,294,299,303,305,313,324,326,327,339,340,344,350,351,354,358,359,361,362,370,372,373,375,377,381,382,385,386,389,392,396,400,403,404,406,410,412,415,416,425,464,467,470,471,474,476,477,478,484,485,486,488,490,493,497,498,499,501,502,503);
my @h1_internal = (20,22,23,24,25,26,33,35,36,41,43,44,49,50,58,59,66,67,69,72,74,75,76,77,78,79,80,81,82,84,87,93,95,96,97,98,104,105,107,108,109,112,117,118,119,121,122,125,128,133,134,140,143,149,151,152,160,161,163,164,165,166,167,174,177,183,189,190,191,192,193,194,195,196,197,198,199,204,208,213,215,216,217,218,219,226,231,233,241,242,243,245,246,247,248,250,254,256,258,259,260,262,263,264,265,266,267,270,271,272,273,280,282,284,293,295,296,297,298,300,301,302,307,308,309,315,316,319,320,323,325,328,330,331,333,334);

my @h3_surface = (25,37,38,41,43,47,48,49,54,61,62,64,66,71,73,79,94,97,98,99,107,108,110,112,120,120,140,142,144,147,148,149,151,153,156,158,159,160,161,173,174,175,176,178,179,181,183,187,188,189,204,205,206,208,209,214,215,224,228,230,238,240,241,255,256,277,278,279,280,285,287,289,292,294,301,305,307,326,328,329,340,341,352,353,356,361,363,364,372,374,375,376,377,383,384,387,391,394,398,402,403,405,406,414,416,418,427,466,466,480,488,491,492,499,500,501,503,505,506,509,510,513,517,518);
my @h3_internal = (27,29,31,32,33,35,42,44,45,50,52,53,58,59,60,67,68,72,75,77,80,82,83,84,85,86,87,88,89,92,95,100,102,103,104,105,106,113,114,115,118,123,124,125,126,127,128,129,131,132,133,134,136,141,143,146,155,163,164,166,167,168,169,170,177,180,182,186,192,193,194,195,196,197,198,199,200,201,202,207,211,216,218,219,220,221,222,229,231,236,244,245,246,247,248,251,253,257,259,260,261,263,265,266,267,268,269,270,272,273,274,281,282,283,284,286,288,291,297,298,299,302,303,304,310,311,318,319,321,322,325,330,332,333,335,336,338,349,350,351,355,358,359,362,367,368,369,373,385,386,389,393,396,411,420,423,425,426,428,429,430,432,434,435,436,437,438,440,441,444,445,446,448,449,452,453,454,455,456,457,458,459,460,463,464,467,469,471,474,476,477,481,483,485,486,487,489,493,494,497,502,508,511,512,515);


my @h1_leading_kr = (4,11,13,16,52,60,73,74,86,88,90,91,97,99,111,113,128,144,151,156,157,162,169,170,171,172,173,178,182,184,199,202,203,205,206,207,209,210,232,240,261,268,269,283,287,289,290,293,326,361,415,488,489);
my @h1_trailing_kr = (3,6,7,11,52,53,64,89,91,99,101,111,129,137,142,144,148,150,156,157,158,162,165,169,172,178,179,185,186,195,197,199,200,201,203,207,231,236,243,251,253,269,274,278,289,290,293,299,324,331,361,389,390,398,415,422,455,467,470,489,510,514,515,526,562);
my @h3_leading_kr = (27,35,57,82,89,94,107,115,138,153,156,163,165,167,169,172,175,176,177,187,188,190,191,192,195,204,218,221,222,224,225,229,234,249,251,254,258,259,293,294,307,308,310,379,393,407,418,482,484,561);
my @h3_trailing_kr = (19,26,27,29,32,59,65,77,79,80,81,82,83,85,88,89,107,117,120,124,126,128,138,144,160,163,169,170,172,174,182,189,191,192,195,196,203,204,205,206,224,225,226,231,233,234,239,241,246,248,252,254,255,257,258,261,276,280,292,301,305,308,311,323,355,358,393,407,417,418,450,458,482,500,521,532,554,562,579);
my @n2_leading_kr = (18,20,23,30,52,93,143,150,194,197,199,208,216,220,221,249,265,307,308,310,313,328,336,339,344,346,368,369,370,372,381,385,387,390,432);
my @n2_trailing_kr = (2,4,5,9,27,30,40,44,45,50,56,65,77,82,83,120,127,147,148,149,151,155,210,216,220,238,248,251,258,263,265,269,302,307,309,310,312,328,329,334,335,338,339,342,347,372,386,392,400,402,403,414,416,432,433,434,455,464);
my @n1_leading_kr = (15,17,23,34,45,64,70,78,105,173,200,214,222,234,248,249,250,254,270,274,275,287,329,332,336,339,344,352,354,367,369,382,390,396,418,427,430,434,451);
my @n1_trailing_kr = (15,17,21,23,38,39,40,42,45,47,48,52,57,67,68,70,73,77,81,82,83,93,100,101,114,130,147,149,188,200,249,254,259,262,264,267,270,273,275,329,331,340,346,352,364,366,367,390,416,418,419,427,435,452,453,455,462);


my @h3_leading_neva = (18,176,159,175,243,242,162,260,264,161,377,213,202,19,292,69,208,178,531,138,172);
my @h3_trailing_neva = (97,79,236,342,538,13,179,468,289,264,235,173,400);
my @n2_leading_neva = (126,56,249,332,399,264,431,215,43,267,248,220,290,328,313,46,208,432,172,372,69,308,370);
my @n2_trailing_neva = (401,372,155,335,220,339,432,430,44,127,263);

my @h1_leading_neva = (4,113,171,257,178,138,290,52,415,184,223,16,102,86,199,157,162,361,200,111,74,145,221,169,64);
my @h1_trailing_neva = (91,200,232,169,268);
my @n1_leading_neva = (163,263,388,6,149,59,78,14,80,101,427,386,200);
my @n1_trailing_neva = (434,275,15,267,83,287);


my @h3_deps_evolving = qw(10 61 151 161 171 174 245 264 347);
## 90 240 277 179 - only in egg-adapted before 1979
my @h1_jianpeng_evolving = qw(156 169 171 203 206 210 238 90 240 277 179);

my @h1_wenfu_evolving = qw(98 110 157 178 202 203 238 176);
#my @not_in_pdb_h1 = (1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 504, 505, 506, 507, 508, 509, 510, 511, 512, 513, 514, 515, 516, 517, 518, 519, 520, 521, 522, 523, 524, 525, 526, 527, 528, 529, 530, 531, 532, 533, 534, 535, 536, 537, 538, 539, 540, 541, 542, 543, 544, 545, 546, 547, 548, 549, 550, 551, 552, 553, 554, 555, 556, 557, 558, 559, 560, 561, 562, 563);
#my @h3_best_sites = (235,202,183,478,538,189,66);
#logic();



#print_tree_with_mutations(176, "n2");
#push @h1_surface, @h1_internal;
#logic_medstat_groups("h1", "nsyn", 1000, "shadow", \@h1_surface);
#logic_medstat_groups("h1", "nsyn", 1000, "leading", \@h1_leading_kr);
#logic_medstat_groups_labelshuffler("h1", "nsyn", 1000, "trailing", \@h1_trailing_kr);
#logic_medstat_groups("h1", "nsyn", 1000, "trailing", \@h1_trailing_kr);

#logic_medstat_groups("h3", "nsyn", 1000, "leading_inter", \@h3_leading_neva);
#logic_medstat_groups_labelshuffler("h3", "nsyn", 1000, "leading_inter", \@h3_leading_neva);
#logic_medstat_groups("h3", "nsyn", 1000, "trailing_inter", \@h3_trailing_neva);
#logic_medstat_groups_labelshuffler("h3", "nsyn", 1000, "trailing_inter", \@h3_trailing_neva);

#logic_medstat_groups("n2", "nsyn", 1000, "leading_inter", \@n2_leading_neva);
#logic_medstat_groups_labelshuffler("n2", "nsyn", 1000, "leading_inter", \@n2_leading_neva);
#logic_medstat_groups("n2", "nsyn", 1000, "trailing_inter", \@n2_trailing_neva);
#logic_medstat_groups_labelshuffler("n2", "nsyn", 1000, "trailing_inter", \@n2_trailing_neva);

#logic_medstat_groups("h1", "nsyn", 1000, "leading_inter", \@h1_leading_neva);
#logic_medstat_groups_labelshuffler("h1", "nsyn", 1000, "leading_inter", \@h1_leading_neva);
#logic_medstat_groups("h1", "nsyn", 1000, "trailing_inter", \@h1_trailing_neva);
#logic_medstat_groups_labelshuffler("h1", "nsyn", 1000, "trailing_inter", \@h1_trailing_neva);

#logic_medstat_groups("n1", "nsyn", 1000, "leading_inter", \@n1_leading_neva);
#logic_medstat_groups_labelshuffler("n1", "nsyn", 1000, "leading_inter", \@n1_leading_neva);
#logic_medstat_groups("n1", "nsyn", 1000, "trailing_inter", \@n1_trailing_neva);
#logic_medstat_groups_labelshuffler("n1", "nsyn", 1000, "trailing_inter", \@n1_trailing_neva);

#logic_medstat_groups("n1", "nsyn", 1000, "leading", \@n1_leading_kr);
#logic_medstat_groups_labelshuffler("n1", "nsyn", 1000, "leading", \@n1_leading_kr);
#logic_medstat_groups("n1", "nsyn", 1000, "trailing", \@n1_trailing_kr);
#logic_medstat_groups_labelshuffler("n1", "nsyn", 1000, "trailing", \@n1_trailing_kr);

#logic_medstat_groups_labelshuffler("n1", "nsyn", 1000, "internal", \@n1_internal);
#logic_medstat_groups("n1", "nsyn", 1000, "internal_dd_norm", \@n1_internal);
#logic_medstat_groups_labelshuffler("n2", "nsyn", 1000, "internal", \@n2_internal);
#logic_medstat_groups("n2", "nsyn", 1000, "internal_dd_norm", \@n2_internal);

#logic_medstat_groups_labelshuffler("h3", "nsyn", 1000, "surface", \@h3_surface);
#logic_medstat_groups("h3", "nsyn", 1000, "surface_dd_norm", \@h3_surface);
#logic_medstat_groups_labelshuffler("h1", "nsyn", 1000, "internal", \@h1_internal);
#logic_medstat_groups("h3", "nsyn", 1000, "internal_dd_norm", \@h3_internal);

#logic_medstat_groups_labelshuffler("n2", "nsyn", 1000, "internal", \@n2_internal);
#logic_medstat_groups("n2", "nsyn", 1000, "internal_dd_norm", \@n2_internal);

#logic_medstat_groups_labelshuffler("n1", "nsyn", 1000, "epitope", \@n1_epitopes);
#logic_medstat_groups_labelshuffler("n2", "nsyn", 1000, "epitope", \@n2_epitopes);
#logic_medstat_groups_labelshuffler("h1", "nsyn", 1000, "epitope", \@h1_epitopes);
#logic_medstat_groups_labelshuffler("h3", "nsyn", 1000, "epitope", \@h3_epitopes);


#logic_medstat_groups_labelshuffler("h1", "nsyn", 1000, "incrbinding", \@h1_increased_binding);
#logic_medstat_groups_labelshuffler("n1", "nsyn", 1000, "wan_epitopes", \@n1_wan_epitopes);
#logic_medstat_groups_labelshuffler("h1", "nsyn", 1000, "host001", \@h1_host_shift_001);
#logic_medstat_groups_labelshuffler("h3", "nsyn", 1000, "host001", \@h3_host_shift_001);
#logic_medstat_groups_labelshuffler("n1", "nsyn", 1000, "host001", \@n1_host_shift_001);
#logic_medstat_groups_labelshuffler("n2", "nsyn", 1000, "host001", \@n2_host_shift_001);

#logic_medstat_groups_labelshuffler("h3", "nsyn", 1000, "smith", \@h3_antigenic_smith);
#logic_medstat_groups_labelshuffler("h3", "nsyn", 1000, "koel", \@h3_antigenic_koel);
#logic_medstat_groups_labelshuffler("n2", "nsyn", 1000, "surface", \@n2_surface);

#logic_medstat_groups("h1", "nsyn", 1000, "antigenic", \@h1_antigenic);
#logic_medstat_groups("h3", "nsyn", 1000, "antigenic", \@h3_antigenic);
##logic_medstat_groups_labelshuffler("h1", "nsyn", 1000, "antigenic", \@h1_antigenic);
#logic_medstat_groups_labelshuffler("h3", "nsyn", 1000, "antigenic", \@h3_antigenic);

#logic_medstat_groups("n1", "nsyn", 1000, "pocket_distance", \@n1_pocket_closest);
#logic_medstat_groups_labelshuffler("n1", "nsyn", 1000, "pocket_distance", \@n1_pocket_closest);

#logic_medstat_groups("h1", "nsyn", 1000, "not_in_pdb", \@h1_pocket_closest);

#logic_medstat_groups("h1", "nsyn", 1000, "epitope_dd_norm", \@h1_epitopes);
#logic_medstat_groups("h3", "nsyn", 1000, "epitope_dd_norm", \@h3_epitopes);
#logic_medstat_groups("n1", "nsyn", 1000, "epitope_dd_norm", \@n1_epitopes);
#logic_medstat_groups("n2", "nsyn", 1000, "epitope_dd_norm", \@n2_epitopes);






#print "h1 nsyn\n";
#logic_global_median_statistics("h3", "nsyn", 0, "koelhist", \@h3_antigenic_koel);
#print "h1 syn\n";
#logic_global_median_statistics("h1", "syn", 0, "dele");
#print "h3 nsyn\n";
#logic_global_median_statistics("h3", "nsyn", 0, "dele");
#print "h3 syn\n";
#logic_global_median_statistics("h3", "syn", 0, "dele");
#print "n1 nsyn\n";
#logic_global_median_statistics("n1", "nsyn", 0, "dele");
#print "n1 syn\n";
#logic_global_median_statistics("n1", "syn", 0, "dele");

##print "n2 nsyn\n";
#logic_global_median_statistics("n2", "nsyn", 0, "dele");
#print "n2 syn\n";
#logic_global_median_statistics("n2", "syn", 0, "dele");




#logic_global_median_statistics("h1", "nsyn", 1000, "norm");
#logic_global_median_statistics("h3", "syn", 1000, "norm");
#logic_global_median_statistics("h1", "syn", 1000, "norm");
#logic_global_median_statistics("h3", "nsyn", 1000, "norm");
#logic_global_median_statistics("n1", "syn", 1000, "norm");
#logic_global_median_statistics("n2", "syn", 1000, "norm");
#logic_global_median_statistics("n1", "nsyn", 1000, "norm");
#logic_global_median_statistics("n2", "nsyn", 1000, "norm");

#logic_median_statistics("h1", "syn", 1000);
#logic_median_statistics("h3", "syn", 1000);
#logic_median_statistics("h1", "nsyn", 1000, "norm");
#logic_median_statistics("h3", "nsyn", 1000);
#logic_median_statistics("n1", "syn", 1000);
#logic_median_statistics("n2", "syn", 1000);
#logic_median_statistics("n1", "nsyn", 1000);
#logic_median_statistics("n2", "nsyn", 1000);
#onemtest(128);
#onemtest(214);
#onemtest(136);
#logic_radius_compare();
#logic_radius_compare_except_seq_and_sis();
#logic_sites_seqsis();
#logic_codon_groups();
#logic_shuffler();
#logic_unrestricted();
#logic();
#logic_collector_general_shuffler("C:/Users/weidewind/Documents/CMD/Coevolution/Influenza/Kryazhimsky11/h3.l.r.newick", "C:/Users/weidewind/Documents/CMD/Coevolution/Influenza/Kryazhimsky11/h3.all.fa");
#logic_collector_sites("C:/Users/weidewind/Documents/CMD/Coevolution/Influenza/Kryazhimsky11/n2.l.r.newick", "C:/Users/weidewind/Documents/CMD/Coevolution/Influenza/Kryazhimsky11/n2.all.fa");
#logic_rcesas_sites("C:/Users/weidewind/Documents/CMD/Coevolution/Influenza/Kryazhimsky11/h3.l.r.newick", "C:/Users/weidewind/Documents/CMD/Coevolution/Influenza/Kryazhimsky11/h3.all.fa");
#print "==============";
#logic("C:/Users/weidewind/Documents/CMD/Coevolution/Influenza/Kryazhimsky11/h3.l.r.newick", "C:/Users/weidewind/Documents/CMD/Coevolution/Influenza/Kryazhimsky11/h3.all.fa");

#logic_general_same_ancestor();
#logic_general_same_ancestor_unrestricted();
sub logic{
my $tree = parse_tree($_[0]);
my %fasta = parse_fasta($_[1]);
my @mutmaps = synmutmap($tree, \%fasta);
my %subs_on_node = %{$mutmaps[0]};
my %nodes_with_sub = %{$mutmaps[1]};


my $wilcox_test = Statistics::Test::WilcoxonRankSum->new();
#470
for (my $ind = 1; $ind <566; $ind++){
	#compute_bitvectors($tree, \%subs_on_node, $ind);
my @distr = find_all_distances($tree, \@{$nodes_with_sub{$ind}}, \%subs_on_node, $ind, 1);

print "Distribution for ".$ind.": ";
print "\nSame aa: ";
	foreach my $dist(@{$distr[0]}){
		print ($dist."\t")
	};
	print "\n";
	print "Different aa: ";
	foreach my $dist(@{$distr[1]}){
		print ($dist."\t")
	};
print "\n";
  try {$wilcox_test->load_data(\@{$distr[0]}, \@{$distr[1]});
  my $prob = $wilcox_test->probability();
  my $pf = sprintf '%f', $prob; # prints 0.091022
  print "\n";
  print $wilcox_test->probability_status();
  print $wilcox_test->summary();
    print "\n";
  }
}



}


sub logic_unrestricted{
my $tree = parse_tree("C:/Users/weidewind/Documents/CMD/Coevolution/Influenza/Kryazhimsky11/h3.l.r.newick");
my %fasta = parse_fasta("C:/Users/weidewind/Documents/CMD/Coevolution/Influenza/Kryazhimsky11/h3.all.fa");

my @ks = keys %fasta;
my $length = length($fasta{$ks[0]});
my @mutmaps = mutmap($tree, \%fasta);
my %subs_on_node = %{$mutmaps[0]};
my %nodes_with_sub = %{$mutmaps[1]};

my @general_same;
my @general_diff;

#my @sieved = sieve(\@{$nodes_with_sub{10}},\%subs_on_node, 10, "Y", 1);
#print "Sieved: ";
#foreach my $s(@sieved){
#	print $$s->get_name()."\t";
#}
my $wilcox_test = Statistics::Test::WilcoxonRankSum->new();
#470
#566
print (($length/3)."\n");
for (my $ind = 1; $ind <=($length/3); $ind++){

#print ("Nodes with sub: ".scalar @{$nodes_with_sub{$ind}});

my @distr = find_min_distances_unrestricted($tree, \@{$nodes_with_sub{$ind}}, \%subs_on_node, $ind, 1);

push @general_same, @{$distr[0]};
push @general_diff, @{$distr[1]};

print "Distribution for ".$ind.": ";
print "\nSame aa: ";
	foreach my $dist(@{$distr[0]}){
		print ($dist."\t")
	};
	print "\n";
	print "Different aa: ";
	foreach my $dist(@{$distr[1]}){
		print ($dist."\t")
	};
print "\n";
  try {$wilcox_test->load_data(\@{$distr[0]}, \@{$distr[1]});
  my $prob = $wilcox_test->probability();
  my $pf = sprintf '%f', $prob; # prints 0.091022
  print "\n";
  print $wilcox_test->probability_status();
  print $wilcox_test->summary();
    print "\n";
  }
  

  
}

print "Same:\n";
foreach my $same(@general_same){
	print $same."\n";
}
print "Different:\n";
foreach my $diff(@general_diff){
	print $diff."\n";
}
  try {$wilcox_test->load_data(\@general_same, \@general_diff);
  my $prob = $wilcox_test->probability();
  my $pf = sprintf '%f', $prob; # prints 0.091022
  print "\n";
  print $wilcox_test->probability_status();
  print $wilcox_test->summary();
    print "\n";
  }

}

## 30.04: for ancestral codon A in site i: (1) normalize summ of distances between the same muts by n of branches that have such mut;
## (2) same normalization - for summ of distances between different muts
## divide 1 by 2
## summ for all sites which have ancestral codon A

sub logic_3 {
		my $tree = parse_tree("C:/Users/weidewind/Documents/CMD/Coevolution/Influenza/Kryazhimsky11/n1.l.r.newick");
		my %fasta = parse_fasta("C:/Users/weidewind/Documents/CMD/Coevolution/Influenza/Kryazhimsky11/n1.all.fa");
		my @mutmaps = synmutmap($tree, \%fasta);
		my %subs_on_node = %{$mutmaps[0]};
		my %nodes_with_sub = %{$mutmaps[1]};

		my @general;

		my $r;
		#print ("MEANTEST: ".mean_ignore_nulls(1,5,$r,6,18));
		my $codon = "CCG";
		
		my $wilcox_test = Statistics::Test::WilcoxonRankSum->new();
	#470
	
	for (my $ind = 1; $ind <566; $ind++){
		
	#compute_bitvectors($tree, \%subs_on_node, $ind);
	#print ("Nodes with sub: ".scalar @{$nodes_with_sub{$ind}});
		my %distr = find_all_distances_3($tree, \@{$nodes_with_sub{$ind}}, \%subs_on_node, $ind, 1);
print "site number $ind";
print "\n";
foreach my $codon(keys %distr){

	my $diff = mean_ignore_nulls(($distr{$codon})->[1]);
	if ($diff == 0) {
		$diff = 0.0001;
	}

	push @general, mean_ignore_nulls(($distr{$codon})->[0])/$diff;
}
	}

foreach my $number(@general){
	print $number."\n";
}

	
}


sub logic_4 {
		my $tree = parse_tree("C:/Users/weidewind/Documents/CMD/Coevolution/Influenza/Kryazhimsky11/n1.l.r.newick");
		my %fasta = parse_fasta("C:/Users/weidewind/Documents/CMD/Coevolution/Influenza/Kryazhimsky11/n1.all.fa");
		my @mutmaps = synmutmap($tree, \%fasta);
		my %subs_on_node = %{$mutmaps[0]};
		my %nodes_with_sub = %{$mutmaps[1]};

		my @general;

		my $r;
		#print ("MEANTEST: ".mean_ignore_nulls(1,5,$r,6,18));
		my $codon = "CCG";
		
		my $wilcox_test = Statistics::Test::WilcoxonRankSum->new();
	#470
	
	for (my $ind = 1; $ind <566; $ind++){
		
	#compute_bitvectors($tree, \%subs_on_node, $ind);
	#print ("Nodes with sub: ".scalar @{$nodes_with_sub{$ind}});
		my %distr = find_all_distances_3($tree, \@{$nodes_with_sub{$ind}}, \%subs_on_node, $ind, 1);
print "site number $ind";
print "\n";
foreach my $codon(keys %distr){

	my $diff = mean_ignore_nulls(($distr{$codon})->[1]);
	if ($diff == 0) {
		$diff = 0.0001;
	}

	push @general, mean_ignore_nulls(($distr{$codon})->[0])/$diff;
}
	}

foreach my $number(@general){
	print $number."\n";
}

	
}



sub logic_radius {
		my $tree = parse_tree("C:/Users/weidewind/Documents/CMD/Coevolution/Influenza/Kryazhimsky11/n1.l.r.newick");
		my %fasta = parse_fasta("C:/Users/weidewind/Documents/CMD/Coevolution/Influenza/Kryazhimsky11/n1.all.fa");
		my @mutmaps = synmutmap($tree, \%fasta);
		my %subs_on_node = %{$mutmaps[0]};
		my %nodes_with_sub = %{$mutmaps[1]};

		my $step = 10;
		my %bins;

		
		my $wilcox_test = Statistics::Test::WilcoxonRankSum->new();
	
	for (my $ind = 1; $ind <566; $ind++){
print "site number $ind";
print "\n";
		my %distr = find_all_distances_radius($tree, \@{$nodes_with_sub{$ind}}, \%subs_on_node, $ind, 1);

foreach my $subst(keys %distr){
	print "subst $subst \n";
	if (defined $distr{$subst}->[0] && defined $distr{$subst}->[1]) {
	my $same_size = $distr{$subst}->[2]->[0];
	my $diff_size = $distr{$subst}->[2]->[1];
	if ($same_size > 1 && $diff_size > 0){
		my %stemp;
		my %dtemp;
		foreach my $s_distance (@{$distr{$subst}->[0]}){
			$stemp{int($s_distance/$step)} = $stemp{int($s_distance/$step)}+1;
			print "Same: distance $s_distance, interval ".int($s_distance/$step)." count ".$stemp{int($s_distance/$step)};
			print "\n";
		}
		foreach my $d_distance (@{$distr{$subst}->[1]}){
			$dtemp{int($d_distance/$step)} = $dtemp{int($d_distance/$step)}+1;
			print "Diff: distance $d_distance, interval ".int($d_distance/$step)." count ".$dtemp{int($d_distance/$step)};
			print "\n";
		}
	
		foreach my $interval((0..50)){
			my $count = $stemp{$interval};
			if (!defined $count){
				$count = 0;
			}
			my $r =  $count/($same_size*($same_size-1));
						print "count: $count, same size: $same_size result $r \n";
			push @{$bins{$interval}->[0]},  $count/($same_size*($same_size-1));  
		}
		foreach my $interval((0..50)){
			my $count =$dtemp{$interval};
						if (!defined $count){
				$count = 0;
			}
			my $r = $count/($same_size*$diff_size);
			print "count: $count, same size: $same_size, diff size: $diff_size result $r \n";
			push @{$bins{$interval}->[1]}, $count/($same_size*$diff_size); 
		}
	}
	}

	
}
	}

foreach my $interval(sort { $a <=> $b } keys %bins){
print "$interval \t";
print mean(@{$bins{$interval}->[0]})."\t";
#print "\n";
#foreach my $number (@{$bins{$interval}->[0]}){
#	print $number."\n";
#}
print mean(@{$bins{$interval}->[1]})."\n";
#print "\n";
#foreach my $number (@{$bins{$interval}->[1]}){
#	print $number."\n";
#}

}
	
}

# for each interval prints two numbers (for each substitution, i.e. same ancestor and same derived aa or codon): observed number of convergent mutations and expected from the of divergent mutations in this interval
# each subst (a, d) is analyzed separately
sub logic_radius_compare {
		my $tree = parse_tree("C:/Users/weidewind/Documents/CMD/Coevolution/Influenza/Kryazhimsky11/n1.l.r.newick");
		my %fasta = parse_fasta("C:/Users/weidewind/Documents/CMD/Coevolution/Influenza/Kryazhimsky11/n1.all.fa");
		my @mutmaps = codonmutmap($tree, \%fasta);
		my %subs_on_node = %{$mutmaps[0]};
		my %nodes_with_sub = %{$mutmaps[1]};

		my $step = 10;
		my %bins;

		
		my $wilcox_test = Statistics::Test::WilcoxonRankSum->new();
	#566
	for (my $ind = 1; $ind <566; $ind++){
print "site number $ind";
print "\n";
		#my %distr = find_all_distances_radius($tree, \@{$nodes_with_sub{$ind}}, \%subs_on_node, $ind, 1);
		my %distr = find_all_distances_radius($tree, \@{$nodes_with_sub{$ind}}, \%subs_on_node, $ind, 1);

foreach my $subst(keys %distr){
	print "subst $subst \n";
	if (defined $distr{$subst}->[0] && defined $distr{$subst}->[1]) {
	my $same_size = $distr{$subst}->[2]->[0];
	my $diff_size = $distr{$subst}->[2]->[1];
	if ($same_size > 0 && $diff_size > 0){
		my %stemp;
		my %dtemp;
		foreach my $s_distance (@{$distr{$subst}->[0]}){
			$stemp{int($s_distance/$step)} = $stemp{int($s_distance/$step)}+1;
			print "Same: distance $s_distance, interval ".int($s_distance/$step)." count ".$stemp{int($s_distance/$step)};
			print "\n";
		}
		foreach my $d_distance (@{$distr{$subst}->[1]}){
			$dtemp{int($d_distance/$step)} = $dtemp{int($d_distance/$step)}+1;
			print "Diff: distance $d_distance, interval ".int($d_distance/$step)." count ".$dtemp{int($d_distance/$step)};
			print "\n";
		}
	
my $obs;
my $exp;	
		foreach my $interval((0..50)){
			my $scount = $stemp{$interval}/$same_size; # normalize to get average count of mutations in i-radius of a mutation
			my $dcount = $dtemp{$interval}/$same_size; 
			if (!defined $scount){
				$scount = 0;
			}
			if (!defined $dcount){
				$dcount = 0;
			}
			$exp +=($same_size-1)*$dcount/$diff_size;
			print "expected ".($same_size-1)*$dcount/$diff_size."\n";
			$obs += $scount;
			print "observed ".$scount."\n";
			push @{$bins{$interval}->[0]}, ($same_size-1)*$dcount/$diff_size;  # expected mean number of convergent mutations in i-radius
			push @{$bins{$interval}->[1]}, $scount; # observed mean number
		}
	if ($exp ne $obs){
		print "Discrepance found! obs $obs, exp $exp \n";
	}	
	else {
		print "OK: obs $obs, exp $exp \n";
	}

	}
	}

	
}
	}

foreach my $interval(sort { $a <=> $b } keys %bins){
	print "interval: $interval expected ".sum(@{$bins{$interval}->[0]})." in ".scalar @{$bins{$interval}->[0]}." observed ".sum(@{$bins{$interval}->[1]})." in ".scalar @{$bins{$interval}->[1]}."\n";
	for (my $num = 0; $num < scalar @{$bins{$interval}->[0]}; $num++){
		print $bins{$interval}->[0]->[$num]."\t".$bins{$interval}->[1]->[$num];
		print "\n";
	}


#  try {$wilcox_test->load_data(\@{$bins{$interval}->[0]}, \@{$bins{$interval}->[1]});
#  my $prob = $wilcox_test->probability();
#  my $pf = sprintf '%f', $prob; # prints 0.091022
#  print "\n";
#  print $wilcox_test->probability_status();
#  print $wilcox_test->summary();
#    print "\n";
#  }
}
	
}


sub logic_rcesas_sites{
	my $tree = parse_tree($_[0]);
	my %fasta = parse_fasta($_[1]);
	my @mutmaps = codonmutmap($tree, \%fasta);
	my %subs_on_node = %{$mutmaps[0]};
	my %nodes_with_sub = %{$mutmaps[1]};


	my $wilcox_test = Statistics::Test::WilcoxonRankSum->new();
	#470
	for (my $ind = 1; $ind <566; $ind++){
		#compute_bitvectors($tree, \%subs_on_node, $ind);
		my %distr = find_all_distances_except_seq_and_sis_radius($tree, \@{$nodes_with_sub{$ind}}, \%subs_on_node, $ind, 1);

	my @same;
	my @diff;
	foreach my $subst(keys %distr){
		print "subst $subst \n";
		if (defined $distr{$subst}->[0] && defined $distr{$subst}->[1]) {
			print " all defined, pushing \n";
			push @same, @{$distr{$subst}->[0]};
			push @diff, @{$distr{$subst}->[1]};
		}
	}

	print "Distribution for ".$ind.": ";
	print "\nSame aa: ";
	foreach my $dist(@same){
		print ($dist."\t")
	};
	print "\n";
	print "Different aa: ";
	foreach my $dist(@diff){
		print ($dist."\t")
	};
	print "\n";
	
  try {$wilcox_test->load_data(\@same, \@diff);
  my $prob = $wilcox_test->probability();
  my $pf = sprintf '%f', $prob; # prints 0.091022
  print "\n";
  print $wilcox_test->probability_status();
  print $wilcox_test->summary();
    print "\n";
  }
}
	
}

sub logic_collector_sites{
	my $tree = parse_tree($_[0]);
	my %fasta = parse_fasta($_[1]);
	my @mutmaps = codonmutmap($tree, \%fasta);
	my %subs_on_node = %{$mutmaps[0]};
	my %nodes_with_sub = %{$mutmaps[1]};


	my $wilcox_test = Statistics::Test::WilcoxonRankSum->new();
	#470
	for (my $ind = 1; $ind <566; $ind++){
		#compute_bitvectors($tree, \%subs_on_node, $ind);
		my %distr = collector($tree, \@{$nodes_with_sub{$ind}}, \%subs_on_node, $ind, 1);

	my @same;
	my @diff;
	foreach my $subst(keys %distr){
		print "subst $subst \n";
		if (defined $distr{$subst}->[0] && defined $distr{$subst}->[1]) {
			print " all defined, pushing \n";
			push @same, @{$distr{$subst}->[0]};
			push @diff, @{$distr{$subst}->[1]};
		}
	}

	print "Distribution for ".$ind.": ";
	print "\nSame aa: ";
	foreach my $dist(@same){
		print ($dist."\t")
	};
	print "\n";
	print "Different aa: ";
	foreach my $dist(@diff){
		print ($dist."\t")
	};
	print "\n";
	
  try {$wilcox_test->load_data(\@same, \@diff);
  my $prob = $wilcox_test->probability();
  my $pf = sprintf '%f', $prob; # prints 0.091022
  print "\n";
  print $wilcox_test->probability_status();
  print $wilcox_test->summary();
    print "\n";
  }
}
	
}

sub logic_collector_general_shuffler {
	my $tree = parse_tree($_[0]);
	my %fasta = parse_fasta($_[1]);
	my @mutmaps = codonmutmap($tree, \%fasta);
	my %subs_on_node = %{$mutmaps[0]};
	my %nodes_with_sub = %{$mutmaps[1]};

	for (my $i = 0; $i <1000; $i++){
	my $ttest = new Statistics::TTest;  
	
	my @same;
	my @diff;
	#470
	for (my $ind = 1; $ind <566; $ind++){
		#compute_bitvectors($tree, \%subs_on_node, $ind);
		my %distr = collector($tree, \@{$nodes_with_sub{$ind}}, \%subs_on_node, $ind, 1);

	foreach my $subst(keys %distr){
		if (defined $distr{$subst}->[0] && defined $distr{$subst}->[1]) {
			push @same, @{$distr{$subst}->[0]};
			push @diff, @{$distr{$subst}->[1]};
		}
	}

}

	
  try {
  	$ttest ->load_data(\@same, \@diff);
    print "$ttest->f_statistic()";
    print "\n";
  }
	}
}

sub logic_collector_general {
	my $tree = parse_tree($_[0]);
	my %fasta = parse_fasta($_[1]);
	my @mutmaps = codonmutmap($tree, \%fasta);
	my %subs_on_node = %{$mutmaps[0]};
	my %nodes_with_sub = %{$mutmaps[1]};


	my $wilcox_test = Statistics::Test::WilcoxonRankSum->new();
	
	my @same;
	my @diff;
	#470
	for (my $ind = 1; $ind <566; $ind++){
		#compute_bitvectors($tree, \%subs_on_node, $ind);
		my %distr = collector($tree, \@{$nodes_with_sub{$ind}}, \%subs_on_node, $ind, 1);

	foreach my $subst(keys %distr){
		print "subst $subst \n";
		if (defined $distr{$subst}->[0] && defined $distr{$subst}->[1]) {
			print " all defined, pushing \n";
			push @same, @{$distr{$subst}->[0]};
			push @diff, @{$distr{$subst}->[1]};
		}
	}

}

	print "\nSame aa: ";
	foreach my $dist(@same){
		print ($dist."\n")
	};
	print "\n";
	print "Different aa: ";
	foreach my $dist(@diff){
		print ($dist."\n")
	};
	print "\n";
	
  try {$wilcox_test->load_data(\@same, \@diff);
  my $prob = $wilcox_test->probability();
  my $pf = sprintf '%f', $prob; # prints 0.091022
  print "\n";
  print $wilcox_test->probability_status();
  print $wilcox_test->summary();
    print "\n";
  }
}


# for each interval prints two numbers (for each substitution, i.e. same ancestor and same derived aa or codon): observed number of convergent mutations and expected from the of divergent mutations in this interval
# each subst (a, d) is analyzed separately
# counts each pair only once
sub logic_radius_compare_except_seq_and_sis {
	my $tree = parse_tree("C:/Users/weidewind/Documents/CMD/Coevolution/Influenza/Kryazhimsky11/h3.l.r.newick");
		my %fasta = parse_fasta("C:/Users/weidewind/Documents/CMD/Coevolution/Influenza/Kryazhimsky11/h3.all.fa");
		my @mutmaps = codonmutmap($tree, \%fasta);
		my %subs_on_node = %{$mutmaps[0]};
		my %nodes_with_sub = %{$mutmaps[1]};
		my $neva = 1;
		my $step = 10;
		my %bins;

		
		my $wilcox_test = Statistics::Test::WilcoxonRankSum->new();
	#566
	for (my $ind = 183; $ind <184; $ind++){
#print "site number $ind";
#print "\n";
		#my %distr = find_all_distances_radius($tree, \@{$nodes_with_sub{$ind}}, \%subs_on_node, $ind, 1);
		my %distr = find_all_distances_except_seq_and_sis_radius($tree, \@{$nodes_with_sub{$ind}}, \%subs_on_node, $ind, 1);
my $testers = 0;
my $testerd = 0;
foreach my $subst(keys %distr){
	#print "subst $subst \n";
	if (defined $distr{$subst}->[0] && defined $distr{$subst}->[1]) {
	my $same_size = $distr{$subst}->[2]->[0];
	my $diff_size = $distr{$subst}->[2]->[1];
	if ($same_size > 0 && $diff_size > 0){
		my %stemp;
		my %dtemp;
		foreach my $s_distance (@{$distr{$subst}->[0]}){
			$stemp{int($s_distance/$step)} = $stemp{int($s_distance/$step)}+1;
			#print "Same: distance $s_distance, interval ".int($s_distance/$step)." count ".$stemp{int($s_distance/$step)};
			##print "\n";
		}
		foreach my $d_distance (@{$distr{$subst}->[1]}){
			$dtemp{int($d_distance/$step)} = $dtemp{int($d_distance/$step)}+1;
			#print "Diff: distance $d_distance, interval ".int($d_distance/$step)." count ".$dtemp{int($d_distance/$step)};
			#print "\n";
		}
	
my $obs;
my $exp;	
		foreach my $interval((0..460)){
			my $scount = $stemp{$interval}/$same_size; # normalize to get average count of mutations in i-radius of a mutation
			my $dcount = $dtemp{$interval}/$same_size; 
			if (!defined $scount){
				$scount = 0;
			}
			if (!defined $dcount){
				$dcount = 0;
			}
			#$exp +=($same_size-1)*$dcount/$diff_size;
			#print "expected ".($same_size-1)*$dcount/$diff_size."\n";
			#$obs += $scount;
			#print "observed ".$scount."\n";
			if ($neva == 1){
				$testerd = $testerd + $dcount/$diff_size;
				$testers = $testers + $scount/($same_size-1);
				push @{$bins{$interval}->[0]}, $dcount/$diff_size;  # expected mean percent of convergent mutations in i-radius
				push @{$bins{$interval}->[1]}, $scount/($same_size-1); # observed mean percent
			}
			else {
				push @{$bins{$interval}->[0]}, ($same_size-1)*$dcount/$diff_size;  # expected mean number of convergent mutations in i-radius
				push @{$bins{$interval}->[1]}, $scount; # observed mean number
			}
		}
	#if ($exp ne $obs){
	#	print "Discrepance found! obs $obs, exp $exp \n";
	#}	
	#else {
	#	print "OK: obs $obs, exp $exp \n";
	#}

	#print "testerd $testerd\n";
	#print "testers $testers\n";
	}
	}

	
}
	}

foreach my $interval(sort { $a <=> $b } keys %bins){
	#print "interval: $interval expected ".sum(@{$bins{$interval}->[0]})." in ".scalar @{$bins{$interval}->[0]}." observed ".sum(@{$bins{$interval}->[1]})." in ".scalar @{$bins{$interval}->[1]}."\n";
	print "$interval,".sum(@{$bins{$interval}->[0]}).",".sum(@{$bins{$interval}->[1]})."\n";
	
	#for (my $num = 0; $num < scalar @{$bins{$interval}->[0]}; $num++){
	#	print $bins{$interval}->[0]->[$num]."\t".$bins{$interval}->[1]->[$num];
	#	print "\n";
	#}


#  try {$wilcox_test->load_data(\@{$bins{$interval}->[0]}, \@{$bins{$interval}->[1]});
#  my $prob = $wilcox_test->probability();
#  my $pf = sprintf '%f', $prob; # prints 0.091022
#  print "\n";
#  print $wilcox_test->probability_status();
#  print $wilcox_test->summary();
#    print "\n";
#  }
}
	
}



# analyze sites
sub logic_median_statistics {
	my $tree = parse_tree("C:/Users/weidewind/Documents/CMD/Coevolution/Influenza/Kryazhimsky11/".$_[0].".l.r.newick");
	my %fasta = parse_fasta("C:/Users/weidewind/Documents/CMD/Coevolution/Influenza/Kryazhimsky11/".$_[0].".all.fa");
	my @mutmaps;
	if($_[1] eq "syn"){
		 @mutmaps = synmutmap($tree, \%fasta);
	} 
	elsif($_[1] eq "nsyn"){
		 @mutmaps = codonmutmap($tree, \%fasta);
	} 
	else {
		die "only syn or nsyn can be used as the second argument; unknown $_[1] was used instead";
	}
	my %subs_on_node = %{$mutmaps[0]};
	my %nodes_with_sub = %{$mutmaps[1]};
	my $step = 1;
	my $iterate = $_[2];
		my $outfile = $_[0]."_".$_[1]."_medstat".$_[3];
	open OUT, ">$outfile" or die "cannot create output file $outfile: $!";
	my @bootstrap_median_diff;
	print OUT "site\tsame_median\tdiff_median\tmedian_difference\tpvalue\n";
#my @arr = (183, 202, 213, 235, 538);
#foreach my $ind(@arr){
for (my $ind = 1; $ind <566; $ind++){
	my @bootstrap_median_diff;
	print OUT $ind."\t";

	my %distr = find_all_distances_except_seq_and_sis_radius($tree, \@{$nodes_with_sub{$ind}}, \%subs_on_node, $ind, 1, $_[1]);
	my @bins = distr_to_stathist(\%distr, $step);
	my $same_median = hist_median(\@{$bins[0]});
	my $diff_median = hist_median(\@{$bins[1]});
	my $obs_difference = $diff_median-$same_median;
	print OUT $same_median."\t"; #this have to be the median of "same" statistics
	print OUT $diff_median."\t"; #this have to be the median of "diff" statistics
	print OUT $obs_difference."\t";
	
	for (my $t = 0; $t < $iterate; $t++){
		my %shuffled_distr = shuffler($tree, \@{$nodes_with_sub{$ind}}, \%subs_on_node, $ind, 1, $_[1]);
		my @shuffler_bins = distr_to_stathist(\%shuffled_distr, $step);
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
	print OUT $pvalue;
	if ($pvalue < 0.01){
		print OUT "\tSignif";
	}
	print OUT "\n";
	}
	close OUT;
}

#global analysis

sub logic_global_median_statistics{
	my $tree = parse_tree("C:/Users/weidewind/Documents/CMD/Coevolution/Influenza/Kryazhimsky11/".$_[0].".l.r.newick");
	my %fasta = parse_fasta("C:/Users/weidewind/Documents/CMD/Coevolution/Influenza/Kryazhimsky11/".$_[0].".all.fa");
	my @mutmaps;
	if($_[1] eq "syn"){
		 @mutmaps = synmutmap($tree, \%fasta);
	} 
	elsif($_[1] eq "nsyn"){
		 @mutmaps = codonmutmap($tree, \%fasta);
	} 
	else {
		die "only syn or nsyn can be used as the second argument; unknown $_[1] was used instead";
	}
	my %subs_on_node = %{$mutmaps[0]};
	my %nodes_with_sub = %{$mutmaps[1]};
	my $iterate = $_[2];
	my $outfile = $_[0]."_".$_[1]."_global_medstat_".$_[3];
	my @array = @{$_[4]};
	open OUT, ">$outfile" or die "cannot create output file $outfile: $!";
	my $step = 10;
	my @bootstrap_median_diff;
	my @bins;
	
if (!@array){
	@array = (1..565);
}
print OUT "site\tsame_median\tdiff_median\tmedian_difference\tpvalue\n";
foreach my $ind(@array){
	#for (my $ind = 1; $ind <566; $ind++){
		print OUT $ind."\n";
		my %distr = find_all_distances_except_seq_and_sis_radius($tree, \@{$nodes_with_sub{$ind}}, \%subs_on_node, $ind, 1, $_[1]);
		my @site_bins = distr_to_stathist(\%distr, $step);

		if (defined $site_bins[0]->[1] && defined $site_bins[1]->[1]){
		for( my $interval = 0; $interval < scalar @{$site_bins[0]}; $interval++){
			$bins[0]->[$interval] += $site_bins[0]->[$interval];
			#print OUT $site_bins[0]->[$interval];
			#print OUT "\t";
			#print OUT $bins[0]->[$interval];
			#print OUT "\t";
		}
		print OUT "\n";
		for (my $interval = 0; $interval < scalar @{$site_bins[1]}; $interval++){
			$bins[1]->[$interval] += $site_bins[1]->[$interval];
			#print OUT $site_bins[1]->[$interval];
			#print OUT "\t";			
			#print OUT $bins[1]->[$interval];
			#print OUT "\t";
		}
		#print OUT "\n";
		}
	}
	
	print OUT "Same\n";
	for( my $interval = 0; $interval < scalar @{$bins[0]}; $interval++){
		print OUT $interval."\t".$bins[0]->[$interval]."\n";
	}
	print OUT "Different\n";
	for( my $interval = 0; $interval < scalar @{$bins[1]}; $interval++){
		print OUT $interval."\t".$bins[1]->[$interval]."\n";
	}
	my $same_median = hist_median(\@{$bins[0]});
	my $diff_median = hist_median(\@{$bins[1]});
	my $obs_difference = $diff_median-$same_median;
	print OUT $same_median."\t"; #this have to be the median of "same" statistics
	print OUT $diff_median."\t"; #this have to be the median of "diff" statistics
	print OUT $obs_difference."\t";
	
	my @bootstrap_median_diff;
	for (my $t = 0; $t < $iterate; $t++){
		my @shuffler_bins;
		for (my $ind = 1; $ind <566; $ind++){
			my %shuffled_distr = shuffler($tree, \@{$nodes_with_sub{$ind}}, \%subs_on_node, $ind, 1, $_[1]);
			my @shuffler_site_bins = distr_to_stathist(\%shuffled_distr, $step);
				if (defined $shuffler_site_bins[0]->[1] && defined $shuffler_site_bins[1]->[1]){

				for(my $interval = 0; $interval < scalar @{$shuffler_site_bins[0]}; $interval++){
					$shuffler_bins[0]->[$interval] += $shuffler_site_bins[0]->[$interval];
				}
				for(my $interval = 0; $interval < scalar @{$shuffler_site_bins[1]}; $interval++){
					$shuffler_bins[1]->[$interval] += $shuffler_site_bins[1]->[$interval];
				}
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
	print OUT $pvalue."\n";
	close OUT;
}


#global analysis

sub global_entrenchment_statistics{
	my $tree = parse_tree("C:/Users/weidewind/Documents/CMD/Coevolution/Influenza/Kryazhimsky11/".$_[0].".l.r.newick");
	my %fasta = parse_fasta("C:/Users/weidewind/Documents/CMD/Coevolution/Influenza/Kryazhimsky11/".$_[0].".all.fa");
	my @mutmaps;
	if($_[1] eq "syn"){
		 @mutmaps = synmutmap($tree, \%fasta);
	} 
	elsif($_[1] eq "nsyn"){
		 @mutmaps = codonmutmap($tree, \%fasta);
	} 
	else {
		die "only syn or nsyn can be used as the second argument; unknown $_[1] was used instead";
	}
	my %subs_on_node = %{$mutmaps[0]};
	my %nodes_with_sub = %{$mutmaps[1]};
	my $iterate = $_[2];
	my $outfile = $_[0]."_".$_[1]."_global_entrenchment_".$_[3];
	my @array = @{$_[4]};
	open OUT, ">$outfile" or die "cannot create output file $outfile: $!";
	my $step = 10;
	my @bootstrap_median_diff;
	my @bins;
	
if (!@array){
	@array = (1..565);
}
print OUT "site\tsame_median\tdiff_median\tmedian_difference\tpvalue\n";
foreach my $ind(@array){
	#for (my $ind = 1; $ind <566; $ind++){
		print OUT $ind."\n";
		my %distr = find_all_distances_except_seq_and_sis_radius($tree, \@{$nodes_with_sub{$ind}}, \%subs_on_node, $ind, 1, $_[1]);
		my @site_bins = distr_to_stathist(\%distr, $step);

		if (defined $site_bins[0]->[1] && defined $site_bins[1]->[1]){
		for( my $interval = 0; $interval < scalar @{$site_bins[0]}; $interval++){
			$bins[0]->[$interval] += $site_bins[0]->[$interval];
			#print OUT $site_bins[0]->[$interval];
			#print OUT "\t";
			#print OUT $bins[0]->[$interval];
			#print OUT "\t";
		}
		print OUT "\n";
		for (my $interval = 0; $interval < scalar @{$site_bins[1]}; $interval++){
			$bins[1]->[$interval] += $site_bins[1]->[$interval];
			#print OUT $site_bins[1]->[$interval];
			#print OUT "\t";			
			#print OUT $bins[1]->[$interval];
			#print OUT "\t";
		}
		#print OUT "\n";
		}
	}
	
	print OUT "Same\n";
	for( my $interval = 0; $interval < scalar @{$bins[0]}; $interval++){
		print OUT $interval."\t".$bins[0]->[$interval]."\n";
	}
	print OUT "Different\n";
	for( my $interval = 0; $interval < scalar @{$bins[1]}; $interval++){
		print OUT $interval."\t".$bins[1]->[$interval]."\n";
	}
	my $same_median = hist_median(\@{$bins[0]});
	my $diff_median = hist_median(\@{$bins[1]});
	my $obs_difference = $diff_median-$same_median;
	print OUT $same_median."\t"; #this have to be the median of "same" statistics
	print OUT $diff_median."\t"; #this have to be the median of "diff" statistics
	print OUT $obs_difference."\t";
	
	my @bootstrap_median_diff;
	for (my $t = 0; $t < $iterate; $t++){
		my @shuffler_bins;
		for (my $ind = 1; $ind <566; $ind++){
			my %shuffled_distr = shuffler($tree, \@{$nodes_with_sub{$ind}}, \%subs_on_node, $ind, 1, $_[1]);
			my @shuffler_site_bins = distr_to_stathist(\%shuffled_distr, $step);
				if (defined $shuffler_site_bins[0]->[1] && defined $shuffler_site_bins[1]->[1]){

				for(my $interval = 0; $interval < scalar @{$shuffler_site_bins[0]}; $interval++){
					$shuffler_bins[0]->[$interval] += $shuffler_site_bins[0]->[$interval];
				}
				for(my $interval = 0; $interval < scalar @{$shuffler_site_bins[1]}; $interval++){
					$shuffler_bins[1]->[$interval] += $shuffler_site_bins[1]->[$interval];
				}
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
	print OUT $pvalue."\n";
	close OUT;
}



## wrong idea, correct implementation - tests medians of distances, not of our statistics.
sub logic_median_test {
		my $tree = parse_tree("C:/Users/weidewind/Documents/CMD/Coevolution/Influenza/Kryazhimsky11/n1.l.r.newick");
		my %fasta = parse_fasta("C:/Users/weidewind/Documents/CMD/Coevolution/Influenza/Kryazhimsky11/n1.all.fa");
		my @mutmaps = synmutmap($tree, \%fasta);
		my %subs_on_node = %{$mutmaps[0]};
		my %nodes_with_sub = %{$mutmaps[1]};
		my $iterate = 1000;
		
print OUT "site\tmedian(same)-median(diff)\tpvalue\n";
#my @arr = (154,158,234,239);
#foreach my $ind(@arr){
	
for (my $ind = 1; $ind <566; $ind++){
	my @same;
	my @diff;
	my @bootstrap_median_diff;

	print "$ind";
	print "\t";

	my %distr = find_all_distances_except_seq_and_sis_radius($tree, \@{$nodes_with_sub{$ind}}, \%subs_on_node, $ind, 1);

	foreach my $subst(keys %distr){

		if (defined $distr{$subst}->[0] && defined $distr{$subst}->[1]) {
			foreach my $s(@{$distr{$subst}->[0]}){
				push @same, $s;

			}
			foreach my $d(@{$distr{$subst}->[1]}){
				push @diff, $d;

			}
		}
	}

	
	for (my $t = 0; $t < $iterate; $t++){
		my @shsame;
		my @shdiff;
		my %shuffled_distr = shuffler($tree, \@{$nodes_with_sub{$ind}}, \%subs_on_node, $ind, 1);
		foreach my $subst(keys %shuffled_distr){
		
			if (defined $shuffled_distr{$subst}->[0] && defined $shuffled_distr{$subst}->[1]) {
				foreach my $s(@{$shuffled_distr{$subst}->[0]}){
					push @shsame, $s;
				}
				foreach my $d(@{$shuffled_distr{$subst}->[1]}){
					push @shdiff, $d;
				}
			}
		}

		#push @bootstrap_median_diff, median_difference(\@shdiff,\@shsame);
		push @bootstrap_median_diff, median(@shdiff)-median(@shsame);
	}
	
		#my $obs_difference = median_difference(\@diff,\@same);
		my $obs_difference = median(@diff)-median(@same);
		print $obs_difference."\t";
		my @sorted = sort {$a <=> $b} @bootstrap_median_diff;
		my $pvalue = 0;
		for (my $i = 0; $i < $iterate; $i++){
			if($sorted[$i] >= $obs_difference){
				$pvalue = ($iterate - $i)/$iterate;
				last;
			}
		}
		print $pvalue;
		#$final{$ind} = $pvalue;
		print "\n";	

	}


#foreach my $site(sort { $a <=> $b } keys %final){
#	print $site."\t".$final{$site}."\n";
#}
		
		
}


sub logic_median_global_test {
		my $tree = parse_tree("C:/Users/weidewind/Documents/CMD/Coevolution/Influenza/Kryazhimsky11/h3.l.r.newick");
		my %fasta = parse_fasta("C:/Users/weidewind/Documents/CMD/Coevolution/Influenza/Kryazhimsky11/h3.all.fa");
		my @mutmaps = codonmutmap($tree, \%fasta);
		my %subs_on_node = %{$mutmaps[0]};
		my %nodes_with_sub = %{$mutmaps[1]};
		my $iterate = 10;
		
	my @same;
	my @diff;
	my @bootstrap_median_diff;
		
for (my $ind = 1; $ind <566; $ind++){



	my %distr = find_all_distances_except_seq_and_sis_radius($tree, \@{$nodes_with_sub{$ind}}, \%subs_on_node, $ind, 1);

	foreach my $subst(keys %distr){

		if (defined $distr{$subst}->[0] && defined $distr{$subst}->[1]) {
			foreach my $s(@{$distr{$subst}->[0]}){
				push @same, $s;
			}
			foreach my $d(@{$distr{$subst}->[1]}){
				push @diff, $d;
			}
		}
	}	
}

	for (my $t = 0; $t < $iterate; $t++){
		my @shsame;
		my @shdiff;
		for (my $ind = 1; $ind <566; $ind++){
			my %shuffled_distr = shuffler($tree, \@{$nodes_with_sub{$ind}}, \%subs_on_node, $ind, 1);
			foreach my $subst(keys %shuffled_distr){
		
				if (defined $shuffled_distr{$subst}->[0] && defined $shuffled_distr{$subst}->[1]) {
					foreach my $s(@{$shuffled_distr{$subst}->[0]}){
						push @shsame, $s;
					}
					foreach my $d(@{$shuffled_distr{$subst}->[1]}){
						push @shdiff, $d;
					}
				}
			}		
		}
		push @bootstrap_median_diff, median(@shdiff)-median(@shsame);
		
	}

		my $obs_difference = median(@diff)-median(@same);
		print $obs_difference."\t";
		my @sorted = sort {$a <=> $b} @bootstrap_median_diff;
		my $pvalue = 0;
		for (my $i = 0; $i < $iterate; $i++){
			if($sorted[$i] >= $obs_difference){
				$pvalue = ($iterate - $i)/$iterate;
				last;
			}
		}
		print $pvalue;
		print "\n";	
		
}


sub test_distr_to_stathist_norm {
	logic_global_median_statistics("h3", "nsyn", 0, "distrtest");
}

sub distr_to_stathist {
	my %distr = %{$_[0]};
	my $step = $_[1];
	my $interval_count = 450/$step;
	my @bins;
	my %hash;
	my %pruned_distr;
	
	# neva july normalization
	foreach my $subst(keys %distr){
		if (defined $distr{$subst}->[0] && defined $distr{$subst}->[1]) {
			my $same_size = $distr{$subst}->[2]->[0];
			my $diff_size = $distr{$subst}->[2]->[1];
			if ($same_size > 0 && $diff_size > 0){
				$pruned_distr{$subst} = $distr{$subst};
				my $ancestor_derived = $subst =~ s/[0-9]//gr; 
				$hash{$ancestor_derived} = 1;
			}
		}
	}
	
	my $mutgroups_count = scalar keys %hash;

	foreach my $subst(keys %pruned_distr){
		
				my $same_size = $pruned_distr{$subst}->[2]->[0];
				my $diff_size = $pruned_distr{$subst}->[2]->[1];
			
					
					my %stemp;
					my %dtemp;
					foreach my $s_distance (@{$pruned_distr{$subst}->[0]}){
						$stemp{int($s_distance/$step)+1} = $stemp{int($s_distance/$step)+1}+1;
					}
					foreach my $d_distance (@{$pruned_distr{$subst}->[1]}){
						$dtemp{int($d_distance/$step)+1} = $dtemp{int($d_distance/$step)+1}+1;
					}
		
					foreach my $interval((0..$interval_count)){
						my $scount = $stemp{$interval}/$same_size; # normalize to get average count of mutations in i-radius of a mutation
						my $dcount = $dtemp{$interval}/$same_size; 
						if (!defined $scount){
							$scount = 0; # can be 0 for mutation on one branch
						}
						if (!defined $dcount){
							$dcount = 0;
						}

							$bins[0]->[$interval] += $scount/(($same_size-1)*$mutgroups_count); # $mutgroups_count - for integral over all intervals to be 1
							$bins[1]->[$interval] += $dcount/($diff_size*$mutgroups_count);
					}

	}
	

	return @bins;
}


sub distr_to_stathist_prev {
	my %distr = %{$_[0]};
	my $step = $_[1];
	my $interval_count = 450/$step;
	my @bins;
	my %hash;
	
	# neva july normalization
	foreach my $subst(keys %distr){
		if (defined $distr{$subst}->[0] && defined $distr{$subst}->[1]) {
			my $same_size = $distr{$subst}->[2]->[0];
			my $diff_size = $distr{$subst}->[2]->[1];
			if ($same_size > 0 && $diff_size > 0){
				my $ancestor_derived = $subst =~ s/[0-9]//gr; 
				$hash{$ancestor_derived} = 1;
			}
		}
	}
	
	my $mutgroups_count = scalar keys %hash;
	
	foreach my $subst(keys %distr){
		
			if (defined $distr{$subst}->[0] && defined $distr{$subst}->[1]) {
				my $same_size = $distr{$subst}->[2]->[0];
				my $diff_size = $distr{$subst}->[2]->[1];
			
				if ($same_size > 0 && $diff_size > 0){
					
					my %stemp;
					my %dtemp;
					foreach my $s_distance (@{$distr{$subst}->[0]}){
						$stemp{int($s_distance/$step)} = $stemp{int($s_distance/$step)}+1;
					}
					foreach my $d_distance (@{$distr{$subst}->[1]}){
						$dtemp{int($d_distance/$step)} = $dtemp{int($d_distance/$step)}+1;
					}
		
					foreach my $interval((0..$interval_count)){
						my $scount = $stemp{$interval}/$same_size; # normalize to get average count of mutations in i-radius of a mutation
						my $dcount = $dtemp{$interval}/$same_size; 
						if (!defined $scount){
							$scount = 0;
						}
						if (!defined $dcount){
							$dcount = 0;
						}

							#push @{$bins[0]->[$interval]}, $dcount/$diff_size;  # expected mean number of convergent mutations in i-radius
							#push @{$bins[1]->[$interval]}, $scount/($same_size-1); # observed mean number
							$bins[0]->[$interval] += $scount/($same_size-1); 
							$bins[1]->[$interval] += $dcount/$diff_size;
					}

				}
			}
	}
	return @bins;
		
}

sub testero {
	my @a1 = (2,3,1,4,5,9,9,9);
	my $m = my_median(\@a1);
	print "mymedian $m ";
	my @a2 = (1,1,1,3,5,10,5);
	my $t = median_difference(\@a1, \@a2);
	print "testero $t";
}

#testero();


sub my_median{
my @values = @{$_[0]};	
my $median;
my $mid = int @values/2;
my @sorted_values = sort @values;
if (@values % 2) {
    $median = $sorted_values[ $mid ];
} else {
    $median = ($sorted_values[$mid-1] + $sorted_values[$mid])/2;
} 
return $median;
}


sub medtest {
	my @array = (5,5,5,5);
	my %hash;
	$hash{1} = 2;
	$hash{2} = 4;
	$hash{3} = 4;
	$hash{4} = 10;
	print hist_median(\@array)."\n";
	print hist_median_for_hash(\%hash);
}

#takes an array of probabilities for 0,1,2...
sub hist_median{
	my @hist = @{$_[0]};
	my $summ = sum (@hist);
	my $head = 0;
	my $interval = 0;
	my $median = 0;
	
	while ($head < $summ/2){
		$head += $hist[$interval];
		$median = $interval;
		$interval++;
	}
	
	if ($head == $summ/2){
		$median += 0.5;
	}
#print_hist(\@hist);
	return $median;
}

sub hist_mean {
	my @hist = @{$_[0]};
	my $summ = sum (@hist);
	my $integer;
	for(my $i = 0; $i <scalar @hist; $i++){
		$integer += $i*$hist[$i];
	}
	return $integer/$summ;
}

#takes a hash of probabilities for 0,1,2...
sub hist_median_for_hash{
	my %prehist =  %{$_[0]};
	my @hist;
	my @sorted_keys = sort {$a <=> $b} keys %prehist;
	for (my $i = 1; $i <= $sorted_keys[-1]; $i++){
		if ($prehist{$i}){
			push @hist, $prehist{$i};
		}
		else {
			push @hist, 0;
		}
	}

	return hist_median(\@hist);
}

sub hist_median_for_hash_arr{
	my %prehist =  %{$_[0]};
	my $number = $_[1];
	my @hist;
	my @sorted_keys = sort {$a <=> $b} keys %prehist;
	for (my $i = 1; $i <= $sorted_keys[-1]; $i++){
		if ($prehist{$i}){
			push @hist, $prehist{$i}[$number];
		}
		else {
			push @hist, 0;
		}
	}

	return hist_median(\@hist);
}

sub hist_mean_for_hash{
	my %prehist =  %{$_[0]};
	my @hist;
	my @sorted_keys = sort {$a <=> $b} keys %prehist;
	for (my $i = 1; $i <= $sorted_keys[-1]; $i++){
		if ($prehist{$i}){
			push @hist, $prehist{$i};
		}
		else {
			push @hist, 0;
		}
	}

	return hist_mean(\@hist);
}

sub hist_mean_for_hash_arr{
	my %prehist =  %{$_[0]};
	my $number = $_[1];
	my @hist;
	my @sorted_keys = sort {$a <=> $b} keys %prehist;
	for (my $i = 1; $i <= $sorted_keys[-1]; $i++){
		if ($prehist{$i}){
			push @hist, $prehist{$i}[$number];
		}
		else {
			push @hist, 0;
		}
	}

	return hist_mean(\@hist);
}


## for hist->interval->site_index
sub hist_median_group {
	my @pre_hist = @{$_[0]};
	my @group = @{$_[1]};
	
	my @hist;
	foreach my $ind(@group){
		for (my $interval = 0; $interval < scalar @pre_hist; $interval++){
			$hist[$interval] += $pre_hist[$interval]->[$ind];
		}
	}
	
	return hist_median(\@hist);
}

sub print_hist {
	my @hist = @{$_[0]};

	my $counter = 0;
	print "\n";
	for (my $interval = 0; $interval < scalar @hist; $interval++){
			print $hist[$interval]."\t";
			$counter+=$hist[$interval];
		}
	print "\n";

}

#test_hist_median_group();

sub test_hist_median_group {
	my @hist = (
				[1,2,0,7,8,9],
				[0,2,1,7,9,8],
				[8,7,9,0,2,1],
	);
	my @group1= (3,4,5);
	my @t1= (24,24,3);
	my @group2= (1,0,2);
	my @t2= (3,3,24);
	my  @group3= (1,0,2, 3,5,4);
	my @t3= (27,27,27);
	print hist_median_group(\@hist, \@group1)."\t" ;
	print hist_median(\@t1)."\t";
	print hist_median_group(\@hist, \@group2)."\t" ;
	print hist_median(\@t2)."\t";
	print hist_median_group(\@hist, \@group3)."\t" ;
	print hist_median(\@t3)."\t";
}


sub median_difference{
	my @a1 = @{$_[0]};
	my @a2 = @{$_[1]};
	my $stat = Statistics::Descriptive::Full->new();
	$stat->add_data(@a1);
	my $median1 = $stat->median();

	$stat = Statistics::Descriptive::Full->new();
	$stat->add_data(@a2);
	my $median2 = $stat->median();
	my $diff = $median1-$median2;
	print "median1 $median1 median2 $median2 diff $diff\n";
	return $diff;
}


sub logic_sites_seqsis {
		my $tree = parse_tree("C:/Users/weidewind/Documents/CMD/Coevolution/Influenza/Kryazhimsky11/h1.l.r.newick");
		my %fasta = parse_fasta("C:/Users/weidewind/Documents/CMD/Coevolution/Influenza/Kryazhimsky11/h1.all.fa");
		my @mutmaps = codonmutmap($tree, \%fasta);
		my %subs_on_node = %{$mutmaps[0]};
		my %nodes_with_sub = %{$mutmaps[1]};
		my $neva = 1;
		my $step = 20;
		my $iterate = 10000;
		my %final;

	#566
	my @arr = (398,238);
	foreach my $ind(@arr){
	#for (my $ind = 1; $ind <566; $ind++){
	my %bins;
	print "site number $ind";
	print "\n";
	my %distr = find_all_distances_except_seq_and_sis_radius($tree, \@{$nodes_with_sub{$ind}}, \%subs_on_node, $ind, 1);
#my $testers = 0;
#my $testerd = 0;
	foreach my $subst(keys %distr){
#print "original key $subst\n";		
	#print "subst $subst \n";
	if (defined $distr{$subst}->[0] && defined $distr{$subst}->[1]) {
	my $same_size = $distr{$subst}->[2]->[0];
	my $diff_size = $distr{$subst}->[2]->[1];
#print " original ss $same_size dd $diff_size\n";
	if ($same_size > 0 && $diff_size > 0){
		my %stemp;
		my %dtemp;
		foreach my $s_distance (@{$distr{$subst}->[0]}){
			$stemp{int($s_distance/$step)} = $stemp{int($s_distance/$step)}+1;
			#print "Same: distance $s_distance, interval ".int($s_distance/$step)." count ".$stemp{int($s_distance/$step)};
			##print "\n";
		}
		foreach my $d_distance (@{$distr{$subst}->[1]}){
			$dtemp{int($d_distance/$step)} = $dtemp{int($d_distance/$step)}+1;
			#print "Diff: distance $d_distance, interval ".int($d_distance/$step)." count ".$dtemp{int($d_distance/$step)};
			#print "\n";
		}
	
my $obs;
my $exp;	
		foreach my $interval((0..230)){
			my $scount = $stemp{$interval}/$same_size; # normalize to get average count of mutations in i-radius of a mutation
			my $dcount = $dtemp{$interval}/$same_size; 
			if (!defined $scount){
				$scount = 0;
			}
			if (!defined $dcount){
				$dcount = 0;
			}
			#$exp +=($same_size-1)*$dcount/$diff_size;
			#print "expected ".($same_size-1)*$dcount/$diff_size."\n";
			#$obs += $scount;
			#print "observed ".$scount."\n";
			if ($neva == 1){
			#	$testerd = $testerd + $dcount/$diff_size;
			#	$testers = $testers + $scount/($same_size-1);
#if ($interval == 0) {
#print "pushed $dcount/$diff_size and $scount/($same_size-1)\n";	
#}
				push @{$bins{$interval}->[0]}, $dcount/$diff_size;  # expected mean number of convergent mutations in i-radius
				push @{$bins{$interval}->[1]}, $scount/($same_size-1); # observed mean number
			}
			else {
				push @{$bins{$interval}->[0]}, ($same_size-1)*$dcount/$diff_size;  # expected mean number of convergent mutations in i-radius
				push @{$bins{$interval}->[1]}, $scount; # observed mean number
			}
		}
	#if ($exp ne $obs){
	#	print "Discrepance found! obs $obs, exp $exp \n";
	#}	
	#else {
	#	print "OK: obs $obs, exp $exp \n";
	#}

	#print "testerd $testerd\n";
	#print "testers $testers\n";
	}
	}

	
}


	my %shuffler_bins; 
	for (my $t = 0; $t < $iterate; $t++){
		my %shuffled_distr = shuffler($tree, \@{$nodes_with_sub{$ind}}, \%subs_on_node, $ind, 1);
		foreach my $subst(keys %shuffled_distr){
#print "shuffled key $subst\n";			
			if (defined $shuffled_distr{$subst}->[0] && defined $shuffled_distr{$subst}->[1]) {
				my $same_size = $shuffled_distr{$subst}->[2]->[0];
				my $diff_size = $shuffled_distr{$subst}->[2]->[1];
#print " shuffled ss $same_size dd $diff_size\n";				
				if ($same_size > 0 && $diff_size > 0){
					my %stemp;
					my %dtemp;
					foreach my $s_distance (@{$shuffled_distr{$subst}->[0]}){
						$stemp{int($s_distance/$step)} = $stemp{int($s_distance/$step)}+1;
					}
					foreach my $d_distance (@{$shuffled_distr{$subst}->[1]}){
						$dtemp{int($d_distance/$step)} = $dtemp{int($d_distance/$step)}+1;
					}
		
					foreach my $interval((0..230)){
						my $scount = $stemp{$interval}/$same_size; # normalize to get average count of mutations in i-radius of a mutation
						my $dcount = $dtemp{$interval}/$same_size; 
						if (!defined $scount){
							$scount = 0;
						}
						if (!defined $dcount){
							$dcount = 0;
						}
						if ($neva == 1){
#if ($interval == 0) {
#print "added $scount/($same_size-1) - $dcount/$diff_size\n";	
#}							
							$shuffler_bins{$interval}->[$t] += $scount/($same_size-1) - $dcount/$diff_size; # same-diff						
						}
						else {
							$shuffler_bins{$interval}->[$t] += $scount - ($same_size-1)*$dcount/$diff_size; # same-diff
						}
					}

				}
			}
		}
	
	}
	
	my $upper;
	foreach my $interval(sort { $a <=> $b } keys %shuffler_bins){
		my $diff = sum(@{$bins{$interval}->[1]})-sum(@{$bins{$interval}->[0]});
		print "interval: $interval \t".$diff."\t";
		my @sorted = sort {$a <=> $b} @{$shuffler_bins{$interval}};
		my $pvalue = 0;
		for (my $i = 0; $i < $iterate; $i++){
			if($sorted[$i] >= $diff){
				$pvalue = ($iterate - $i)/$iterate;
				last;
			}
		}
		print $pvalue;
		$final{$ind}->{$interval} = $pvalue;
		#my $stat = Statistics::Descriptive::Full->new();
		#$stat->add_data(\@{$shuffler_bins{$interval}});
		#my $upper = $stat->percentile(95);
		#my $uppermore = $stat->percentile(99);
		#my $uppermost = $stat->percentile(99.9);
		#my $lower = $stat->percentile(5);
		#print "upper\t$upper\tuppermore\t$uppermore\tuppermost\t$uppermost";
		#print "\t";
		#print "lower\t$lower";
		#if ($upper < $diff) {
		#	print "\tSignificant?";
		#}
		print "\n";	
		
	}




	}


foreach my $site(sort { $a <=> $b } keys %final){
	print $site."\t";
	foreach my $bin (sort { $a <=> $b } keys %{$final{$site}}){
		print $final{$site}->{$bin};
		print "\t";
	}
	print "\n";
}

	
}






#
# 
sub logic_shuffler {
	my $tree = parse_tree("C:/Users/weidewind/Documents/CMD/Coevolution/Influenza/Kryazhimsky11/h3.l.r.newick");
		my %fasta = parse_fasta("C:/Users/weidewind/Documents/CMD/Coevolution/Influenza/Kryazhimsky11/h3.all.fa");
		my @mutmaps = synmutmap($tree, \%fasta);
		my %subs_on_node = %{$mutmaps[0]};
		my %nodes_with_sub = %{$mutmaps[1]};
		my $neva = 1;
		my $step = 1;
		my %bins;

	#566
	for (my $t = 0; $t < 10000; $t++){
	for (my $ind = 1; $ind <566; $ind++){
print "site number $ind";
print "\n";
		#my %distr = find_all_distances_radius($tree, \@{$nodes_with_sub{$ind}}, \%subs_on_node, $ind, 1);
		my %distr = shuffler($tree, \@{$nodes_with_sub{$ind}}, \%subs_on_node, $ind, 1);
my $testers = 0;
my $testerd = 0;
foreach my $subst(keys %distr){
	if (defined $distr{$subst}->[0] && defined $distr{$subst}->[1]) {
	my $same_size = $distr{$subst}->[2]->[0];
	my $diff_size = $distr{$subst}->[2]->[1];
	if ($same_size > 0 && $diff_size > 0){
		my %stemp;
		my %dtemp;
		foreach my $s_distance (@{$distr{$subst}->[0]}){
			$stemp{int($s_distance/$step)} = $stemp{int($s_distance/$step)}+1;
		}
		foreach my $d_distance (@{$distr{$subst}->[1]}){
			$dtemp{int($d_distance/$step)} = $dtemp{int($d_distance/$step)}+1;
		}
		
		foreach my $interval((0..460)){
			my $scount = $stemp{$interval}/$same_size; # normalize to get average count of mutations in i-radius of a mutation
			my $dcount = $dtemp{$interval}/$same_size; 
			if (!defined $scount){
				$scount = 0;
			}
			if (!defined $dcount){
				$dcount = 0;
			}

			if ($neva == 1){
				$testerd = $testerd + $dcount/$diff_size;
				$testers = $testers + $scount/($same_size-1);
				$bins{$interval}->[$t] += $scount/($same_size-1) - $dcount/$diff_size; # same-diff
			}
			else {
				$bins{$interval}->[$t] += $scount - ($same_size-1)*$dcount/$diff_size; # same-diff
			}
		}

	}
	}

	
}
	}
}
foreach my $interval(sort { $a <=> $b } keys %bins){
	print "interval: $interval \n";
	for (my $num = 0; $num < 10000; $num++){
		print $bins{$interval}->[$num];
		print "\n";
	}

}
	
}


# for each interval prints two numbers (for each substitution, i.e. same ancestor and same derived aa or codon): observed number of convergent mutations and expected from the of divergent mutations in this interval
# each subst (a, d) is analyzed separately
# counts each pair only once
sub logic_codon_groups {
	my $tree = parse_tree("C:/Users/weidewind/Documents/CMD/Coevolution/Influenza/Kryazhimsky11/n2.l.r.newick");
		my %fasta = parse_fasta("C:/Users/weidewind/Documents/CMD/Coevolution/Influenza/Kryazhimsky11/n2.all.fa");
		my @mutmaps = synmutmap($tree, \%fasta);
		my %subs_on_node = %{$mutmaps[0]};
		my %nodes_with_sub = %{$mutmaps[1]};
		my %codon_evolution = codon_evolution("n2");

		my $step = 10;
		my %bins_neutral;
		my %bins_changing;
		my $neva = 1;

		
		my $wilcox_test = Statistics::Test::WilcoxonRankSum->new();
	#566
	for (my $ind = 1; $ind <566; $ind++){
print "site number $ind";
print "\n";
		#my %distr = find_all_distances_radius($tree, \@{$nodes_with_sub{$ind}}, \%subs_on_node, $ind, 1);
		my %distr = find_all_distances_codon_groups($tree, \@{$nodes_with_sub{$ind}}, \%subs_on_node, $ind, \%codon_evolution);

foreach my $subst(keys %distr){
	print "subst $subst \n";
	if (defined $distr{$subst}->[0] && defined $distr{$subst}->[1]) {
	my $same_size = $distr{$subst}->[2]->[0];
	my $diff_size = $distr{$subst}->[2]->[1];
	my $different_usefulness = $distr{$subst}->[2]->[2];
	
	if ($same_size > 0 && $diff_size > 0){
		my %stemp;
		my %dtemp;
		foreach my $s_distance (@{$distr{$subst}->[0]}){
			$stemp{int($s_distance/$step)} = $stemp{int($s_distance/$step)}+1;
			print "Same: distance $s_distance, interval ".int($s_distance/$step)." count ".$stemp{int($s_distance/$step)};
			print "\n";
		}
		foreach my $d_distance (@{$distr{$subst}->[1]}){
			$dtemp{int($d_distance/$step)} = $dtemp{int($d_distance/$step)}+1;
			print "Diff: distance $d_distance, interval ".int($d_distance/$step)." count ".$dtemp{int($d_distance/$step)};
			print "\n";
		}
	
my $obs;
my $exp;	
		foreach my $interval((0..50)){
			my $scount = $stemp{$interval}/$same_size; # normalize to get average count of mutations in i-radius of a mutation
			my $dcount = $dtemp{$interval}/$same_size; 
			if (!defined $scount){
				$scount = 0;
			}
			if (!defined $dcount){
				$dcount = 0;
			}
			$exp +=($same_size-1)*$dcount/$diff_size;
			print "expected ".($same_size-1)*$dcount/$diff_size."\n";
			$obs += $scount;
			print "observed ".$scount."\n";
			if ($neva == 1){
				if ($different_usefulness == 0){
					push @{$bins_neutral{$interval}->[0]}, $dcount/$diff_size;  # expected mean number of convergent mutations in i-radius
					push @{$bins_neutral{$interval}->[1]}, $scount/($same_size-1); # observed mean number
				}
				else{
					push @{$bins_changing{$interval}->[0]}, $dcount/$diff_size;  # expected mean number of convergent mutations in i-radius
					push @{$bins_changing{$interval}->[1]}, $scount/($same_size-1); # observed mean number
				}
			}
			else {
				if ($different_usefulness == 0){
					push @{$bins_neutral{$interval}->[0]}, ($same_size-1)*$dcount/$diff_size;  # expected mean number of convergent mutations in i-radius
					push @{$bins_neutral{$interval}->[1]}, $scount; # observed mean number
				}
				else{
					push @{$bins_changing{$interval}->[0]}, ($same_size-1)*$dcount/$diff_size;  # expected mean number of convergent mutations in i-radius
					push @{$bins_changing{$interval}->[1]}, $scount; # observed mean number
				}
			}
		}
	if ($exp ne $obs){
		print "Discrepance found! obs $obs, exp $exp \n";
	}	
	else {
		print "OK: obs $obs, exp $exp \n";
	}

	}
	}

	
}
	}

foreach my $interval(sort { $a <=> $b } keys %bins_neutral){
	print "interval: $interval expected ".sum(@{$bins_neutral{$interval}->[0]})." in ".scalar @{$bins_neutral{$interval}->[0]}." observed ".sum(@{$bins_neutral{$interval}->[1]})." in ".scalar @{$bins_neutral{$interval}->[1]}."\n";
	for (my $num = 0; $num < scalar @{$bins_neutral{$interval}->[0]}; $num++){
		print $bins_neutral{$interval}->[0]->[$num]."\t".$bins_neutral{$interval}->[1]->[$num];
		print "\n";
	}
}

foreach my $interval(sort { $a <=> $b } keys %bins_changing){
	print "Changing\n";
	print "interval: $interval expected ".sum(@{$bins_changing{$interval}->[0]})." in ".scalar @{$bins_changing{$interval}->[0]}." observed ".sum(@{$bins_changing{$interval}->[1]})." in ".scalar @{$bins_changing{$interval}->[1]}."\n";
	for (my $num = 0; $num < scalar @{$bins_changing{$interval}->[0]}; $num++){
		print $bins_changing{$interval}->[0]->[$num]."\t".$bins_changing{$interval}->[1]->[$num];
		print "\n";
	}
}
	
}

# for each interval prints two numbers (for each substitution, i.e. same ancestor and same derived aa or codon): observed number of convergent mutations and expected from the of divergent mutations in this interval
# data for different substs from one site is summed up
sub logic_radius_compare_heap {
		my $tree = parse_tree("C:/Users/weidewind/Documents/CMD/Coevolution/Influenza/Kryazhimsky11/n1.l.r.newick");
		my %fasta = parse_fasta("C:/Users/weidewind/Documents/CMD/Coevolution/Influenza/Kryazhimsky11/n1.all.fa");
		my @mutmaps = codonmutmap($tree, \%fasta);
		my %subs_on_node = %{$mutmaps[0]};
		my %nodes_with_sub = %{$mutmaps[1]};

		my $step = 10;
		my %bins;

		
		my $wilcox_test = Statistics::Test::WilcoxonRankSum->new();
	
	for (my $ind = 1; $ind <566; $ind++){
print "site number $ind";
print "\n";
		my %distr = find_all_distances_radius($tree, \@{$nodes_with_sub{$ind}}, \%subs_on_node, $ind, 1);

my $same_size;
my $diff_size;

		my %stemp;
		my %dtemp;

foreach my $subst(keys %distr){
	print "subst $subst \n";
	if (defined $distr{$subst}->[0] && defined $distr{$subst}->[1]) {
		$same_size += $distr{$subst}->[2]->[0];
		$diff_size += $distr{$subst}->[2]->[1];


		foreach my $s_distance (@{$distr{$subst}->[0]}){
			$stemp{int($s_distance/$step)} = $stemp{int($s_distance/$step)}+1;
			print "Same: distance $s_distance, interval ".int($s_distance/$step)." count ".$stemp{int($s_distance/$step)};
			print "\n";
		}
		foreach my $d_distance (@{$distr{$subst}->[1]}){
			$dtemp{int($d_distance/$step)} = $dtemp{int($d_distance/$step)}+1;
			print "Diff: distance $d_distance, interval ".int($d_distance/$step)." count ".$dtemp{int($d_distance/$step)};
			print "\n";
		}


	
	}

	
}
		if ($same_size >0 && $diff_size>0){
			foreach my $interval((0..50)){
				my $scount = $stemp{$interval}/$same_size; # normalize to get average count of mutations in i-radius of a mutation
				my $dcount = $dtemp{$interval}/$same_size; 
				if (!defined $scount){
					$scount = 0;
				}
				if (!defined $dcount){
					$dcount = 0;
				}
				print $scount/($same_size-1)."\t";
				push @{$bins{$interval}->[0]}, ($same_size-1)*$dcount/$diff_size;  # expected mean number of convergent mutations in i-radius
				push @{$bins{$interval}->[1]}, $scount; # observed mean number
			}
		}

	}

foreach my $interval(sort { $a <=> $b } keys %bins){
	print "interval: $interval \n";
	foreach (my $num = 0; $num < scalar @{$bins{$interval}->[0]}; $num++){
		print $bins{$interval}->[0]->[$num]."\t".$bins{$interval}->[1]->[$num];
		print "\n";
	}

  try {$wilcox_test->load_data(\@{$bins{$interval}->[0]}, \@{$bins{$interval}->[1]});
#  my $prob = $wilcox_test->probability();
#  my $pf = sprintf '%f', $prob; # prints 0.091022
#  print "\n";
#  print $wilcox_test->probability_status();
#  print $wilcox_test->summary();
#    print "\n";
  }
}
	
}

#todo

sub logic_gold {
		my $tree = parse_tree("C:/Users/weidewind/Documents/CMD/Coevolution/Influenza/Kryazhimsky11/n1.l.r.newick");
		my %fasta = parse_fasta("C:/Users/weidewind/Documents/CMD/Coevolution/Influenza/Kryazhimsky11/n1.all.fa");
		my @mutmaps = synmutmap($tree, \%fasta);
		my %subs_on_node = %{$mutmaps[0]};
		my %nodes_with_sub = %{$mutmaps[1]};

		my @general;

		my $r;
		#print ("MEANTEST: ".mean_ignore_nulls(1,5,$r,6,18));
		my $codon = "CCG";
		
		my $wilcox_test = Statistics::Test::WilcoxonRankSum->new();
	#470
	
	for (my $ind = 1; $ind <566; $ind++){

		my %distr = find_all_distances_gold($tree, \@{$nodes_with_sub{$ind}}, \%subs_on_node, $ind, 1);
		print "site number $ind";
		print "\n";

	}
my %pairhash = pairhash();
foreach my $key(%pairhash){
	print "conv: ";
	print $pairhash{$key} -> {"conv_count"};
	print " div: ";
	print $pairhash{$key} -> {"div_count"};
	print " dist: ";
	print $pairhash{$key} -> {"distance"};
	print "\n";
}

	
}

## bootstrap: standard shuffler, shuffles labels on sites
sub logic_medstat_groups_labelshuffler {
	my $tree = parse_tree("C:/Users/weidewind/Documents/CMD/Coevolution/Influenza/Kryazhimsky11/".$_[0].".l.r.newick");
	my %fasta = parse_fasta("C:/Users/weidewind/Documents/CMD/Coevolution/Influenza/Kryazhimsky11/".$_[0].".all.fa");
	my @mutmaps;
	if($_[1] eq "syn"){
		 @mutmaps = synmutmap($tree, \%fasta);
	} 
	elsif($_[1] eq "nsyn"){
		 @mutmaps = codonmutmap($tree, \%fasta);
	} 
	else {
		die "only syn or nsyn can be used as the second argument; unknown $_[1] was used instead";
	}
	my %subs_on_node = %{$mutmaps[0]};
	my %nodes_with_sub = %{$mutmaps[1]};
	my $iterate = $_[2];
		my $outfile = $_[0]."_".$_[1]."_groups_labelshuffler_".$_[3];
	my @group = @{$_[4]};
	my @complement;
	
	open OUT, ">$outfile" or die "cannot create output file $outfile: $!";
	my $step = 1;
	my @bootstrap_median_diff;
	my @bins;
	my @meaningful_sites;
	
	
	
#my @arr = @h1_host_shift;
print OUT "site\tsame_median\tdiff_median\tmedian_difference\tpvalue\n";
#foreach my $ind(@arr){
	for (my $ind = 1; $ind <566; $ind++){
		print OUT $ind."\n";
		my %distr = find_all_distances_except_seq_and_sis_radius($tree, \@{$nodes_with_sub{$ind}}, \%subs_on_node, $ind, 1, $_[1]);
		my @site_bins = distr_to_stathist(\%distr, $step);

		if (defined $site_bins[0]->[1] && defined $site_bins[1]->[1]){
		push @meaningful_sites, $ind;	
		for( my $interval = 0; $interval < scalar @{$site_bins[0]}; $interval++){
			$bins[0]->[$interval]->[$ind] = $site_bins[0]->[$interval];
		}

		for (my $interval = 0; $interval < scalar @{$site_bins[1]}; $interval++){
			$bins[1]->[$interval]->[$ind] = $site_bins[1]->[$interval];
		}

		#}
		}
	}
	
	my %group_hash;
	foreach my $gs(@group){
		$group_hash{$gs} = 1;
	}
	my @group;
	
	foreach my $ms(@meaningful_sites){
		if (exists $group_hash{$ms}){
			push @group, $ms;
		}
		else {
			push @complement, $ms;
		}
	}

	my $same_median_group = hist_median_group(\@{$bins[0]}, \@group);
	my $diff_median_group = hist_median_group(\@{$bins[1]}, \@group);
	my $obs_difference_group = $diff_median_group-$same_median_group;
	print OUT "\t size \t same median \t diff median \t difference\n";
	print OUT "group\t";
	print OUT scalar @group;
	print OUT "\t";
	print OUT $same_median_group."\t"; #this have to be the median of "same" statistics
	print OUT $diff_median_group."\t"; #this have to be the median of "diff" statistics
	print OUT $obs_difference_group."\n";
	
	my $same_median_complement = hist_median_group(\@{$bins[0]}, \@complement);
	my $diff_median_complement = hist_median_group(\@{$bins[1]}, \@complement);
	my $obs_difference_complement = $diff_median_complement-$same_median_complement;
	print OUT "complement\t";
	print OUT scalar @complement;
	print OUT "\t";
	print OUT $same_median_complement."\t"; #this have to be the median of "same" statistics
	print OUT $diff_median_complement."\t"; #this have to be the median of "diff" statistics
	print OUT $obs_difference_complement."\n";
	
	my $diffdiff = $obs_difference_group - $obs_difference_complement;
	print OUT $diffdiff."\n";
	
	my @bootstrap_median_diff;
	my @group_bootstrap;
	
	my $group_count;
	my $enrich_count;
	my $depl_count;
	
	for (my $t = 0; $t < $iterate; $t++){
		my @shuffler_bins;
		
		for (my $ind = 1; $ind <566; $ind++){
	
			my %distr = shuffler($tree, \@{$nodes_with_sub{$ind}}, \%subs_on_node, $ind, 1, $_[1]);
			my @site_bins = distr_to_stathist(\%distr, $step);
			
			if (defined $site_bins[0]->[1] && defined $site_bins[1]->[1]){
			push @meaningful_sites, $ind;	
			for( my $interval = 0; $interval < scalar @{$site_bins[0]}; $interval++){
				$shuffler_bins[0]->[$interval]->[$ind] = $site_bins[0]->[$interval];
			}

			for (my $interval = 0; $interval < scalar @{$site_bins[1]}; $interval++){
				$shuffler_bins[1]->[$interval]->[$ind] = $site_bins[1]->[$interval];
			}

			}
			
		}
		my $bootstrap_difference_group = hist_median_group(\@{$shuffler_bins[1]}, \@group)-hist_median_group(\@{$shuffler_bins[0]}, \@group);
		my $bootstrap_difference_complement = hist_median_group(\@{$shuffler_bins[1]}, \@complement)-hist_median_group(\@{$shuffler_bins[0]}, \@complement);	
		if ($bootstrap_difference_group >= $obs_difference_group){
			$group_count++;
		}
		if ($bootstrap_difference_group-$bootstrap_difference_complement >= $diffdiff){
			$enrich_count++;
		}
		if ($bootstrap_difference_group-$bootstrap_difference_complement <= $diffdiff){
			$depl_count++;
		}

	}
	
	print OUT "group pvalue ".$group_count/$iterate."\n";
	print OUT "enrichment pvalue ".$enrich_count/$iterate."\n";
	print OUT "depletion pvalue ".$depl_count/$iterate."\n";
	close OUT;
}

## bootstrap: randomly chooses group of sites
sub logic_medstat_groups {
	my $tree = parse_tree("C:/Users/weidewind/Documents/CMD/Coevolution/Influenza/Kryazhimsky11/".$_[0].".l.r.newick");
	my %fasta = parse_fasta("C:/Users/weidewind/Documents/CMD/Coevolution/Influenza/Kryazhimsky11/".$_[0].".all.fa");
	my @mutmaps;
	if($_[1] eq "syn"){
		 @mutmaps = synmutmap($tree, \%fasta);
	} 
	elsif($_[1] eq "nsyn"){
		 @mutmaps = codonmutmap($tree, \%fasta);
	} 
	else {
		die "only syn or nsyn can be used as the second argument; unknown $_[1] was used instead";
	}
	my %subs_on_node = %{$mutmaps[0]};
	my %nodes_with_sub = %{$mutmaps[1]};
	my $iterate = $_[2];
	my $outfile = $_[0]."_".$_[1]."_groups_medstat_".$_[3];
	my @group = @{$_[4]};
	my @complement;
	
	
	open OUT, ">$outfile" or die "cannot create output file $outfile: $!";
	my $step = 1;
	my @bootstrap_median_diff;
	my @bins;
	my %sites_hash;
	my @meaningful_sites;
	
#my @arr = (154,158,234,239);
print OUT "site\tsame_median\tdiff_median\tmedian_difference\tpvalue\n";
#foreach my $ind(@arr){
	for (my $ind = 1; $ind <566; $ind++){
		print OUT $ind."\n";
		my %distr = find_all_distances_except_seq_and_sis_radius($tree, \@{$nodes_with_sub{$ind}}, \%subs_on_node, $ind, 1, $_[1]);
		my @site_bins = distr_to_stathist(\%distr, $step);
		if (defined $site_bins[0]->[1] && defined $site_bins[1]->[1]){
		push @meaningful_sites, $ind;	
		for( my $interval = 0; $interval < scalar @{$site_bins[0]}; $interval++){
			$bins[0]->[$interval]->[$ind] = $site_bins[0]->[$interval];
		}

		for (my $interval = 0; $interval < scalar @{$site_bins[1]}; $interval++){
			$bins[1]->[$interval]->[$ind] = $site_bins[1]->[$interval];
		}

		
		}
	}
	
	my %group_hash;
	foreach my $gs(@group){
		$group_hash{$gs} = 1;
	}
	my @group;
	
	foreach my $ms(@meaningful_sites){
		if (exists $group_hash{$ms}){
			push @group, $ms;
		}
		else {
			push @complement, $ms;
		}
	}
	
	
	my $same_median_group = hist_median_group(\@{$bins[0]}, \@group);
	my $diff_median_group = hist_median_group(\@{$bins[1]}, \@group);
	my $obs_difference_group = $diff_median_group-$same_median_group;
	print OUT "\t same median \t diff median \t difference\n";
	print OUT "group\t";
	print OUT $same_median_group."\t"; #this have to be the median of "same" statistics
	print OUT $diff_median_group."\t"; #this have to be the median of "diff" statistics
	print OUT $obs_difference_group."\n";
	
	my $same_median_complement = hist_median_group(\@{$bins[0]}, \@complement);
	my $diff_median_complement = hist_median_group(\@{$bins[1]}, \@complement);
	my $obs_difference_complement = $diff_median_complement-$same_median_complement;
	print OUT "complement\t";
	print OUT $same_median_complement."\t"; #this have to be the median of "same" statistics
	print OUT $diff_median_complement."\t"; #this have to be the median of "diff" statistics
	print OUT $obs_difference_complement."\n";
	
	my $diffdiff = $obs_difference_group - $obs_difference_complement;
	print OUT $diffdiff."\n";
	
	#my @bootstrap_median_diff_group;
	#my @bootstrap_median_diff_complement;
	my $counter1 = 0;
	my $counter2 = 0;
	my $counter3 = 0;
	my $counter4 = 0;
	my $counter5;
	my $counter6;


	for (my $t = 0; $t < $iterate; $t++){
		my @bootstrap_group = shuffle @meaningful_sites;
		my @bootstrap_complement = splice (@bootstrap_group, scalar @group, scalar @meaningful_sites - scalar @group);
		
		my $same_median_group = hist_median_group(\@{$bins[0]}, \@bootstrap_group);
		my $diff_median_group = hist_median_group(\@{$bins[1]}, \@bootstrap_group);
	#print OUT $same_median_group."\t"; #this have to be the median of "same" statistics
	#print OUT $diff_median_group."\t"; #this have to be the median of "diff" statistics	
	#print OUT $diff_median_group-$same_median_group."\n";
		#push @bootstrap_median_diff_group,  $diff_median_group-$same_median_group;
	
		my $same_median_complement = hist_median_group(\@{$bins[0]}, \@bootstrap_complement);
		my $diff_median_complement = hist_median_group(\@{$bins[1]}, \@bootstrap_complement);
	#print OUT $same_median_complement."\t"; #this have to be the median of "same" statistics
	#print OUT $diff_median_complement."\t"; #this have to be the median of "diff" statistics
	#print OUT $diff_median_complement-$same_median_complement."\n";
	#print OUT "___\n";
		#push @bootstrap_median_diff_complement, $diff_median_complement-$same_median_complement;
		if ($diff_median_group-$same_median_group - $diff_median_complement+$same_median_complement >= $diffdiff){
			$counter5++;
		}
		if ($diff_median_group-$same_median_group - $diff_median_complement+$same_median_complement <= $diffdiff){
			$counter6++;
		}
		if ($diff_median_group-$same_median_group >= $obs_difference_group){ 
			$counter1++;
			if ($diff_median_complement-$same_median_complement <= $obs_difference_complement){
				$counter2++;
			}
		}
		
		if ($diff_median_group-$same_median_group <= $obs_difference_group){ 
			$counter3++;
			if ($diff_median_complement-$same_median_complement >= $obs_difference_complement){
				$counter4++;
			}
		}
		
	
	}
	
	print OUT "pvalue e  ".$counter1/$iterate." pvalue enrichment  ".$counter2/$iterate."\n"; 
	print OUT "pvalue d ".$counter3/$iterate." pvalue depletion ".$counter4/$iterate."\n";
	print OUT "pvalue diffdiff enrichment ".$counter5/$iterate." pvalue diffdiff depletion ".$counter6/$iterate."\n";
	close OUT;
}



sub test_group_bootstrap{
		my @meaningful_sites = (1,2,3,6,7,8,11,12,13);
		my @group = (1,2,3);
		my @bootstrap_group = shuffle @meaningful_sites;
		my @bootstrap_complement = splice (@bootstrap_group, scalar @group, scalar @meaningful_sites - scalar @group);
		
		foreach my $g(@bootstrap_group){
			print $g.", ";
		}
		print "\n";
		foreach my $c(@bootstrap_complement){
			print $c.", ";
		}
}


sub logic_general_groups {
my $tree = parse_tree("C:/Users/weidewind/Documents/CMD/Coevolution/Influenza/Kryazhimsky11/n1.l.r.newick");
my %fasta = parse_fasta("C:/Users/weidewind/Documents/CMD/Coevolution/Influenza/Kryazhimsky11/n1.all.fa");
my @mutmaps = mutmap($tree, \%fasta);
my %subs_on_node = %{$mutmaps[0]};
my %nodes_with_sub = %{$mutmaps[1]};
my @ks = keys %fasta;
my $length = length($fasta{$ks[0]});
my @indexes = (1..$length/3);

my $wilcox_test = Statistics::Test::WilcoxonRankSum->new();
#470 n1
#566 h1

#h1 epitopes (checked)
#my @group1 = qw(141 142 171 173 175 176 178 179 180 169 172 205 206 209 211 153 156 158 182 186 195 220 237 238 253 286 87 88 90 91 92 132);

#n1 epitopes
my @group1 = qw(380 381 382 383 384 385 386 388 389 390 393 397 398 199 200 201 202 223 329 330 332 331 333 336 337 339 340 341 343 344 356 363 364 365 366 367);

#h3 epitopes
# letter + 16
#my @group1 = qw( 138 140 142 147 148 146 149 151 153 154 156 158 159 160 161 162 166 168 184 144 145 171 172 173 174 175 176 179 180 181 202 203 204 205 206 208 209 210 212 213 214 60 61 62 63 64 66 67 69 70 289 291 292 294 295 296 310 313 315 316 320 321 323 324 325 326 327 328 112 118 119 133 137 183 186 187 188 189 190 191 192 193 195 198 217 219 223 224 225 228 229 230 231 232 233 234 235 242 243 244 245 246 254 256 258 260 262 263 264 73 75 78 79 83 91 94 96 97 98 99 102 103 104 107 108 110 125 276 277 278 281 );

#n2 epitopes
#my @group1 = qw(383 384 385 386 387 389 390 391 392 393 394 396 399 400 401 403 197 198 199 200 221 222 328 329 330 331 332 334 336 338 339 341 342 343 344 346 347 357 358 359  366 367 368 369 370);

#h1 Huang (antigenic), as is in file Tables_main (Huang + 17, from msa)
#my @group1 = qw( 138 144 145 147 150 158 163 142 170 177 200 203 206 207 208 210 211 52 53 60 288 290 291 294 312 327 111 180 222 226 233 239 241 64 71 86 88 90 97 99 284 );

#n2 surface (as is in the surface txt file)
#my @group1 = qw( 88 89 90 91 92 93 468 469 470 463 464 465 466 455 456 457 450 451 452 453 413 414 415 416 417 399 400 401 402 403 383 384 385 386 387 388 389 390 391 392 393 394 366 367 368 369 370 371 356 357 358 359 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50 51 52 53 54 55 56 57 58 59 60 61 62 63 64 65 66 67 68 69 70 71 72 73 74 75 76 77 78 79 80 81 82 83 84 85 86 341 342 343 344 336 337 338 339 328 329 330 331 332 306 307 308 309 310 311 312 313 283 284 285 286 267 268 269 270 271 261 262 263 264 265 244 245 246 247 248 249 250 251 218 219 220 221 222 208 209 210 196 197 198 199 200 161 162 149 150 151 152 153 154 125 126 127 128 110 111 112 113 461 459 437 396 381 380 378 347 346 334 326 315 304 296 295 292 277 273 259 258 253 236 234 224 216 215 212 189 187 173 171 169 147 146 143 141 130 118 107 95 );


#h3 antigentic Steinbruck
#my @group1 = qw( 138 160 171 223 161 205 233 294 66 153 174 276 140 151 230 278 78 172 212 292 41 91 99 147 202 218 238 241 19 204 69 180 190 209 217 229 246 98 149 159 162 176 213 18 70 188 260 206 242 );
my $group_size = scalar @group1;



my %hash1;
for my $ind(@group1){
	$hash1{$ind} = 1;
}

my %hash_of_dist;

	my @general_same1;
	my @general_diff1;
	my @general_same2;
	my @general_diff2;
	for (my $ind = 1; $ind <=($length/3); $ind++){
		compute_bitvectors($tree, \%subs_on_node, $ind);
		my @distr = find_min_distances_naive($tree, \@{$nodes_with_sub{$ind}}, \%subs_on_node, $ind, 1);
		$hash_of_dist{$ind} = \@distr;
		if (exists $hash1{$ind}){
			push @general_same1, @{$distr[0]};
			push @general_diff1, @{$distr[1]};
		}
		else {
			push @general_same2, @{$distr[0]};
			push @general_diff2, @{$distr[1]};	
		}
	}

print_info_for_hist( \@general_same1, \@general_diff1 );
print_info_for_hist( \@general_same2, \@general_diff2 );

	my @s1 = grep { defined() and length() } @general_same1;
	my @d1 = grep { defined() and length() } @general_diff1;
	my @s2 = grep { defined() and length() } @general_same2;
	my @d2 = grep { defined() and length() } @general_diff2;
	print "\n Checking real data ";
	print median(@s1)."\t".median(@d1)."\t".median(@s2)."\t".median(@d2);
	print "\n";

#push @general_same1, @general_same2;
#push @general_diff1, @general_diff2;
#for (my $i = 0; $i < 50; $i++){
#	my @sample = (shuffle(@general_same1))[0..$group_size-1];
#}

print "Median difference distribution:\n";
for (my $i = 0; $i < 200; $i++){
	my @general_same1;
	
	
	my @general_diff1;
	my @general_same2;
	my @general_diff2;
	my @shuffled = shuffle(@indexes);
	my @sample = @shuffled[0..$group_size-1];
	my @complement = @shuffled[$group_size..$length/3-1];
	
	for my $s(@sample){
		my @distr = @{$hash_of_dist{$s}};
		push @general_same1, @{$distr[0]};
		push @general_diff1, @{$distr[1]};
		
	}
		for my $c(@complement){
		my @distr = @{$hash_of_dist{$c}};
		push @general_same2, @{$distr[0]};
		push @general_diff2, @{$distr[1]};
	}
	my $check0 = scalar @sample;
	my $check00 = scalar @complement;
	my $check1 = scalar @general_same1;
	my $check2 = scalar @general_same2;
	my $check3 = scalar @general_diff1;
	my $check4 = scalar @general_diff2;
	print "sizes: $check0 $check00 $check1 $check2 $check3 $check4";
	print "\n";
	my @s1 = grep { defined() and length() } @general_same1;
	my @d1 = grep { defined() and length() } @general_diff1;
	my @s2 = grep { defined() and length() } @general_same2;
	my @d2 = grep { defined() and length() } @general_diff2;
	print median(@s1)."\t".median(@d1)."\t".median(@s2)."\t".median(@d2);
	print "\n";
}


}


sub print_info_for_hist{
	my @distr1 = @{$_[0]};
	my @distr2 = @{$_[1]};
	my $wilcox_test = Statistics::Test::WilcoxonRankSum->new();
	
	print "Same:\n";
	foreach my $same(@distr1){
		print $same."\n";
	}
	print "Different:\n";
	foreach my $diff(@distr2){
		print $diff."\n";
	}
  	try { 
  		$wilcox_test->load_data(\@distr1, \@distr2);
 	 	my $prob = $wilcox_test->probability();
  		my $pf = sprintf '%f', $prob; # prints 0.091022
  		print "\n";
 		print $wilcox_test->probability_status();
  		print $wilcox_test->summary();
    	print "\n";
  	}
	
}

sub logic_general {
	my $tree = parse_tree("C:/Users/weidewind/Documents/CMD/Coevolution/Influenza/Kryazhimsky11/n1.l.r.newick");
my %fasta = parse_fasta("C:/Users/weidewind/Documents/CMD/Coevolution/Influenza/Kryazhimsky11/n1.all.fa");
my @mutmaps = synmutmap($tree, \%fasta);
my %subs_on_node = %{$mutmaps[0]};
my %nodes_with_sub = %{$mutmaps[1]};

my @general_same;
my @general_diff;

#my @sieved = sieve(\@{$nodes_with_sub{10}},\%subs_on_node, 10, "Y", 1);
#print "Sieved: ";
#foreach my $s(@sieved){
#	print $$s->get_name()."\t";
#}
my $wilcox_test = Statistics::Test::WilcoxonRankSum->new();
#470
for (my $ind = 1; $ind <566; $ind++){
	#compute_bitvectors($tree, \%subs_on_node, $ind);
#print ("Nodes with sub: ".scalar @{$nodes_with_sub{$ind}});
#my @distr = find_min_distances_naive($tree, \@{$nodes_with_sub{$ind}}, \%subs_on_node, $ind, 1);
my @distr = find_all_distances($tree, \@{$nodes_with_sub{$ind}}, \%subs_on_node, $ind, 1);
push @general_same, @{$distr[0]};
push @general_diff, @{$distr[1]};

}
print "Same:\n";
foreach my $same(@general_same){
	print $same."\n";
}
print "Different:\n";
foreach my $diff(@general_diff){
	print $diff."\n";
}
  try {$wilcox_test->load_data(\@general_same, \@general_diff);
  my $prob = $wilcox_test->probability();
  my $pf = sprintf '%f', $prob; # prints 0.091022
  print "\n";
  print $wilcox_test->probability_status();
  print $wilcox_test->summary();
    print "\n";
  }
	
}


## compare overall distribution of min distances between the same mutations and between different mutations in one site.
## Difference from logic_general: only considers mutations with the same ancestor aa (or codon) 
## (_unrestricted) -> ignores oaks between maples

sub logic_general_same_ancestor_unrestricted {
	my $tree = parse_tree("C:/Users/weidewind/Documents/CMD/Coevolution/Influenza/Kryazhimsky11/n1.l.r.newick");
my %fasta = parse_fasta("C:/Users/weidewind/Documents/CMD/Coevolution/Influenza/Kryazhimsky11/n1.all.fa");
my @mutmaps = codonmutmap($tree, \%fasta);
my %subs_on_node = %{$mutmaps[0]};
my %nodes_with_sub = %{$mutmaps[1]};

my @general_same;
my @general_diff;

my $wilcox_test = Statistics::Test::WilcoxonRankSum->new();
#470 n1 468 n2
#566 h1
for (my $ind = 1; $ind <566; $ind++){
	compute_bitvectors($tree, \%subs_on_node, $ind);
#print ("Nodes with sub: ".scalar @{$nodes_with_sub{$ind}});
my @distr = find_min_distances_same_ancestor_unrestricted($tree, \@{$nodes_with_sub{$ind}}, \%subs_on_node, $ind, 1);

push @general_same, @{$distr[0]};
push @general_diff, @{$distr[1]};

}
print "Same:\n";
foreach my $same(@general_same){
	print $same."\n";
}
print "Different:\n";
foreach my $diff(@general_diff){
	print $diff."\n";
}
  try {$wilcox_test->load_data(\@general_same, \@general_diff);
  my $prob = $wilcox_test->probability();
  my $pf = sprintf '%f', $prob; # prints 0.091022
  print "\n";
  print $wilcox_test->probability_status();
  print $wilcox_test->summary();
    print "\n";
  }
	
}


#logic_2();
## mutmap or synmutmap inside - for nsyn or syn substitutions respectively
##control distribution for site i - min distances to any mutation at each site (set size ~ length of the protein*number of branches with mutation in site i)
sub logic_2{
	my $tree = parse_tree("C:/Users/weidewind/Documents/CMD/Coevolution/Influenza/Kryazhimsky11/n1.l.r.newick");
	my %fasta = parse_fasta("C:/Users/weidewind/Documents/CMD/Coevolution/Influenza/Kryazhimsky11/n1.all.fa");
	my @mutmaps = synmutmap($tree, \%fasta);
	my %subs_on_node = %{$mutmaps[0]};
	my %nodes_with_sub = %{$mutmaps[1]};


	my $wilcox_test = Statistics::Test::WilcoxonRankSum->new();

	compute_min_distances_in_subtree($tree, \%subs_on_node);
	my %min_dists_to_any; ## node name to a hash : key = site_index, value = min_dist to any mutation in this site

 	my $root = $tree->get_root();
 	$root->visit_depth_first(
				-post => sub{ #all daughters have been processed
					my $node=shift;
					my %min_dists = compute_min_distances_global($node, \%subs_on_node);
					$min_dists_to_any{$node->get_name} = \%min_dists;
				}
	);

	#470
	for (my $ind = 1; $ind <566; $ind++){
		compute_bitvectors($tree, \%subs_on_node, $ind);
		#print ("Nodes with sub: ".scalar @{$nodes_with_sub{$ind}});
		my @distr = find_min_distances_naive($tree, \@{$nodes_with_sub{$ind}}, \%subs_on_node, $ind, 1);


		print "Distribution for ".$ind.": ";
		print "\nSame aa: ";
		foreach my $dist(@{$distr[0]}){
			print ($dist."\t");
		}
		print "\n";
		
		my @distr_control;
		my @nodes = @{$nodes_with_sub{$ind}};
		for my $n(@nodes){
			push @distr_control, values $min_dists_to_any{$$n->get_name()};
		}
		
		print "Any mutation in any site: ";
		foreach my $dist(@distr_control){
			print ($dist."\t");
		}
		print "\n";
  		try {
  			$wilcox_test->load_data(\@{$distr[0]}, \@distr_control);
 			my $prob = $wilcox_test->probability();
  			my $pf = sprintf '%f', $prob; # prints 0.091022
  			print "\n";
  			print $wilcox_test->probability_status();
  			print $wilcox_test->summary();
   			print "\n";
  		}
	}
}

sub subtreetest{
	my $tree = parse_tree("C:/Users/weidewind/Documents/CMD/Coevolution/Influenza/Kryazhimsky11/h1.l.r.newick");
my %fasta = parse_fasta("C:/Users/weidewind/Documents/CMD/Coevolution/Influenza/Kryazhimsky11/h1.all.fa");
my @mutmaps = mutmap($tree, \%fasta);
my %subs_on_node = %{$mutmaps[0]};
my %nodes_with_sub = %{$mutmaps[1]};

 compute_all_distances_in_subtree($tree, \%subs_on_node);
 my $root = $tree->get_root();
	$root->visit_depth_first(
				-post => sub{ #all daughters have been processed
					my $node=shift;
					my %dists = compute_all_distances_global($node, \%subs_on_node);
					print "node name ";
					print $node->get_name;
					print "\n";
					for my $site(keys %dists){
						print ("site name ".$site);
						print "\n";
						for my $d(@{$dists{$site}}){
							print ($d."_");
							print "\t";
						}
					}
					print "\n";
				}
	);
}


## For every node of the given tree sets a hash -min_distances_in_subtree:
## key = index of site (in protein sequence)
## value = minimal distance to any mutation in this site that happened in the subtree down of this node.
## (to find whole-tree minimal distance, you need such arrays for all nodes on the path from this node to the root)


sub compute_min_distances_global{
	my $node = $_[0];
	my %subs_on_node = %{$_[1]};
	#my $seq_length = $_[2];
	
	my %min_dists = %{$node->get_generic("-min_distances_in_subtree")};

	my $tnode = $node;
	
	#print ("My node name ".$node->get_name()."\t");
	#my @ancestors = $node->get_ancestors();
	#for my $anc(@ancestors){
		while(!$tnode->is_root()){
			my $sister = ${$tnode->get_sisters()}[1]; #some strange manipulations to get the correct () sister (#supposing that every node (except the root) has one and only one sister )
			if ($sister->get_name() eq $tnode->get_name()){
				$sister = ${$tnode->get_sisters()}[0];
			}
		#print ("sister ".$sister->get_name()."\t");
			my %sister_min_dists = %{$sister->get_generic("-min_distances_in_subtree")}; 
		#print (" number of subs ".(scalar keys %sister_min_dists)."\t");	
			my %sister_subs = %{$subs_on_node{$sister->get_name()}};
			my $dist_to_sister = $tnode->get_branch_length() + $sister->get_branch_length();
			for my $site_index(keys %sister_subs){ # if there is a mutation in this site in sister node, check if it's closer than the closest mutation in your own subtree 
				if (!exists $min_dists{$site_index} || $dist_to_sister < $min_dists{$site_index}){
					$min_dists{$site_index} = $dist_to_sister;
				}
			}
			for my $site_index(keys %sister_min_dists ){
				if (!exists $sister_subs{$site_index}){ # if there is no such mutation in sister node, try min distance from the sister subtree
					if(!exists $min_dists{$site_index} || $sister_min_dists{$site_index} + $dist_to_sister < $min_dists{$site_index}){
						$min_dists{$site_index} = $sister_min_dists{$site_index} + $dist_to_sister;
					}
				}
			}
		
			$tnode = $tnode->get_parent();
		}
	#print "\n";
	return %min_dists;
	
}

sub compute_all_distances_global{
	my $node = $_[0];
	my %subs_on_node = %{$_[1]};
	#my $seq_length = $_[2];
	
	my %dists = %{$node->get_generic("-all_distances_in_subtree")};

	my $tnode = $node;
	
	#print ("My node name ".$node->get_name()."\t");
	#my @ancestors = $node->get_ancestors();
	#for my $anc(@ancestors){
		while(!$tnode->is_root()){
			my $sister = ${$tnode->get_sisters()}[1]; #some strange manipulations to get the correct () sister (#supposing that every node (except the root) has one and only one sister )
			if ($sister->get_name() eq $tnode->get_name()){
				$sister = ${$tnode->get_sisters()}[0];
			}
		#print ("sister ".$sister->get_name()."\t");
			my %sister_all_dists = %{$sister->get_generic("-all_distances_in_subtree")}; 
		#print (" number of subs ".(scalar keys %sister_min_dists)."\t");	
			my %sister_subs = %{$subs_on_node{$sister->get_name()}};
			my $dist_to_sister = $tnode->get_branch_length() + $sister->get_branch_length();
			for my $site_index(keys %sister_subs){ # if there is a mutation in this site in sister node, check if it's closer than the closest mutation in your own subtree 
				if (!exists $dists{$site_index}){
								my @dist_array;
								push @dist_array, $dist_to_sister;
								$dists{$site_index} = \@dist_array;
				}
				else {
								my @dist_array = @{$dists{$site_index}};
								push @dist_array, $dist_to_sister;
								$dists{$site_index} = \@dist_array;
								#todo: check 
				}
			}
			for my $site_index(keys %sister_all_dists ){
	 #  distances from the sister subtree
					if(!exists $dists{$site_index}){
						my @dist_array;
						for my $sisdist(@{$sister_all_dists{$site_index}}){
							push @dist_array, $sisdist + $dist_to_sister;
						}
						$dists{$site_index} = \@dist_array;
					}
					else {
								my @dist_array = @{$dists{$site_index}};
								for my $sisdist(@{$sister_all_dists{$site_index}}){
									push @dist_array, $sisdist + $dist_to_sister;
								}
								$dists{$site_index} = \@dist_array;
								#todo: check 
					}
				
			}
		
			$tnode = $tnode->get_parent();
		}
	#print "\n";
	return %dists;
	
}





sub compute_min_distances_in_subtree{
	my $tree = $_[0];
	my %subs_on_node = %{$_[1]};
	
	my $root = $tree->get_root();
	$root->visit_depth_first(
				-post => sub{ #all daughters have been processed
					my $node=shift;
					my $i = 0;
					my %min_dists;
					#print ("Node name: ".$node->get_name()."\t");
					while (my $child = $node->get_child($i)){
						$i++;
						#print ("Child $i name: ".$child->get_name()."\t");
						my %child_min_dists = %{$child->get_generic("-min_distances_in_subtree")};
						foreach my $site(keys %child_min_dists){ # add child branch length to child distances
							my $new_dist = $child_min_dists{$site}+$child->get_branch_length();
							if (!exists $min_dists{$site} || $min_dists{$site} > $new_dist){ # do not change min_dist, if existing distance (from the other child, obviously) is less
								$min_dists{$site} = $new_dist;
							}
						}
						foreach my $new_site(keys %{$subs_on_node{$child->get_name()}}){ # add distances to mutations in the child itself, if they sre less than already computed ones.
							my $new_dist = $child->get_branch_length();
							if (!exists $min_dists{$new_site} || $min_dists{$new_site} > $new_dist){
								$min_dists{$new_site} = $new_dist;
							}
						}
					}
					#print (" number of sites with muutations in the subtree ".(scalar keys %min_dists)."\n");
					$node->set_generic("-min_distances_in_subtree" => \%min_dists);
					#print $node->get_name()."\t";
					#print Dumper (%min_dists);
					#print "\n";
					
				}
			);
}


sub compute_all_distances_in_subtree {
	my $tree = $_[0];
	my %subs_on_node = %{$_[1]};
	
	my $root = $tree->get_root();
	$root->visit_depth_first(
				-post => sub{ #all daughters have been processed
					my $node=shift;
					my $i = 0;
					my %dists;
					print ("Node name: ".$node->get_name()."\t");
					while (my $child = $node->get_child($i)){
						$i++;
						print ("Child $i name: ".$child->get_name()."\t");
						my %child_dists = %{$child->get_generic("-all_distances_in_subtree")};
						foreach my $site(keys %child_dists){ # add child branch length to child distances
							foreach my $new_dist(@{$child_dists{$site}}){
								$new_dist += $child->get_branch_length();
								if (!exists $dists{$site}){
									my @dist_array;
									push @dist_array, $new_dist;
									print " pushing in new $new_dist ";
									$dists{$site} = \@dist_array;
								}
								else {
									my @dist_array = @{$dists{$site}};
									push @dist_array, $new_dist;
									print " pushing in existing $new_dist ";
									$dists{$site} = \@dist_array;
									#todo: check 
								}
							}
						}
						foreach my $new_site(keys %{$subs_on_node{$child->get_name()}}){ # add distances to mutations in the child itself.
							my $new_dist = $child->get_branch_length();
							if (!exists $dists{$new_site}){ 
								my @dist_array;
								push @dist_array, $new_dist;
								print " pushing in new child distance $new_dist ";
								$dists{$new_site} = \@dist_array;
							}
							else {
								my @dist_array = @{$dists{$new_site}};
								push @dist_array, $new_dist;
								print " pushing in existing child distance $new_dist ";
								$dists{$new_site} = \@dist_array;
								#todo: check 
							}
						}
					}
					print (" number of sites with muutations in the subtree ".(scalar keys %dists)."\n");
					$node->set_generic("-all_distances_in_subtree" => \%dists);
					print $node->get_name()."\t";
					#print Dumper (%min_dists);
					#print "\n";
					
				}
			);
}




sub find_min_distances_unrestricted {
	my $tree = $_[0];
	my @nodes = @{$_[1]};
	my %subs_on_node = %{$_[2]};
	my $site_index = $_[3];
	#my $same_derived = $_[4];
	my @min_distances_same;
	my @min_distances_diff;
	
	for (my $i = 0; $i < scalar @nodes; $i++){

		my $derived1 = ${$subs_on_node{${$nodes[$i]}->get_name()}}{$site_index}->{"Substitution::derived_allele"};

		for (my $j = 0; $j < scalar @nodes; $j++){
			if ($j == $i){ next; }
			my $dist = calc_true_patristic_distance(${$nodes[$i]}, ${$nodes[$j]});
			my $derived2 = ${$subs_on_node{${$nodes[$j]}->get_name()}}{$site_index}->{"Substitution::derived_allele"};
			if ($derived1 eq $derived2){
					push @min_distances_same, $dist;
			}
			else {
					push @min_distances_diff, $dist;
				}
		}

	}
	return (\@min_distances_same, \@min_distances_diff);	
}


## collects all distances,  ignores oaks between maples

sub find_all_distances_same_ancestor_unrestricted {
	my $tree = $_[0];
	my @nodes = @{$_[1]};
	my %subs_on_node = %{$_[2]};
	my $site_index = $_[3];
	#my $same_derived = $_[4];
		my $myCodonTable   = Bio::Tools::CodonTable->new();
	my @distances_same;
	my @distances_diff;
	
	for (my $i = 0; $i < scalar @nodes; $i++){

#		my $min_dist_same;
#		my $min_dist_diff;
		
		my $sub1 = ${$subs_on_node{${$nodes[$i]}->get_name()}}{$site_index};
		my $derived1 = $myCodonTable->translate($sub1->{"Substitution::derived_allele"});
		my $ancestor1 = $sub1->{"Substitution::ancestral_allele"};

		for (my $j = 0; $j < scalar @nodes; $j++){
			if ($j == $i){ next; }
			my $sub2 = ${$subs_on_node{${$nodes[$j]}->get_name()}}{$site_index};
			my $derived2 = $myCodonTable->translate($sub2->{"Substitution::derived_allele"});
			my $ancestor2 = $sub2->{"Substitution::ancestral_allele"};
			if ($ancestor1 ne $ancestor2 ){ next; }
			my $dist = calc_true_patristic_distance(${$nodes[$i]}, ${$nodes[$j]});
			if ($derived1 eq $derived2){
				push @distances_same, $dist;
			}
			else {
				push @distances_diff, $dist;
			}
		}

	}
	return (\@distances_same, \@distances_diff);	
}


## collects all distances,  ignores oaks between maples. Same ancestor. Returns arrays. For synmutmap, logic_3

sub find_all_distances {
	my $tree = $_[0];
	my @nodes = @{$_[1]};
	my %subs_on_node = %{$_[2]};
	my $site_index = $_[3];

	my @distances_same;
	my @distances_diff;
	
	for (my $i = 0; $i < scalar @nodes; $i++){
		
		my $sub1 = ${$subs_on_node{${$nodes[$i]}->get_name()}}{$site_index};
		my $derived1 = $sub1->{"Substitution::derived_allele"};
		my $ancestor1 = $sub1->{"Substitution::ancestral_allele"};
			
		for (my $j = 0; $j < scalar @nodes; $j++){
			if ($j == $i){ next; }
			my $sub2 = ${$subs_on_node{${$nodes[$j]}->get_name()}}{$site_index};
			my $derived2 = $sub2->{"Substitution::derived_allele"};
			my $ancestor2 = $sub2->{"Substitution::ancestral_allele"};
			if ($ancestor1 ne $ancestor2 ){ next; }
			my $dist = calc_true_patristic_distance(${$nodes[$i]}, ${$nodes[$j]});

			if ($derived1 eq $derived2){
				push @distances_same, $dist;
			}
			else {
				push @distances_diff, $dist;
			}
		}
	}
	return  (\@distances_same, \@distances_diff);	
}

## 
sub find_all_distances_radius {
	my $tree = $_[0];
	my @nodes = @{$_[1]};
	my %subs_on_node = %{$_[2]};
	my $site_index = $_[3];

	my %hash;
	my @distances_same;
	my @distances_diff;
	my $myCodonTable   = Bio::Tools::CodonTable->new();
	
	for (my $i = 0; $i < scalar @nodes; $i++){

		
		my $sub1 = ${$subs_on_node{${$nodes[$i]}->get_name()}}{$site_index};
		my $derived1 = $myCodonTable -> translate($sub1->{"Substitution::derived_allele"});
	#	my $derived1 = ($sub1->{"Substitution::derived_allele"});
		my $ancestor1 = $sub1->{"Substitution::ancestral_allele"};
			
		my $count_same = 1; # to add the node1 itself
		my $count_diff;
print "node 1:  $ancestor1 $derived1 \n";
		for (my $j = 0; $j < scalar @nodes; $j++){
			if ($j == $i){ next; }
			my $sub2 = ${$subs_on_node{${$nodes[$j]}->get_name()}}{$site_index};
			my $derived2 = $myCodonTable -> translate($sub2->{"Substitution::derived_allele"});
		#	my $derived2 = ($sub2->{"Substitution::derived_allele"});
			my $ancestor2 = $sub2->{"Substitution::ancestral_allele"};
			if ($ancestor1 ne $ancestor2 ){ next; }
			print "node 2:  $ancestor2 $derived2 \n";
			my $dist = calc_true_patristic_distance(${$nodes[$i]}, ${$nodes[$j]});
			if (!exists $hash{"$ancestor1$derived1"} ){
				my @same = ();
				my @diff = ();
				$hash{"$ancestor1$derived1"} = (\@same, \@diff);
			}

			if ($derived1 eq $derived2){
				push @{ ($hash{"$ancestor1$derived1"})->[0] }, $dist;
				$count_same++;
			}
			else {
				push @{ ($hash{"$ancestor1$derived1"})->[1] }, $dist;
				$count_diff++;
			}


		}
		
		print " count same: $count_same, count diff: $count_diff \n";
		push @{ ($hash{"$ancestor1$derived1"})->[2] }, $count_same;
		push @{ ($hash{"$ancestor1$derived1"})->[2] }, $count_diff;
		

	}
	return %hash;	
}




# ignore sequential nodes and immediate sisters; all counts are kept separately for each node
sub find_all_distances_codon_groups {
	my $tree = $_[0];
	my @nodes = @{$_[1]};
	my %subs_on_node = %{$_[2]};
	my $site_index = $_[3];
	my %codon_evolution = %{$_[4]};

	my %hash;
	my @distances_same;
	my @distances_diff;
	my $myCodonTable   = Bio::Tools::CodonTable->new();
	
	for (my $i = 0; $i < scalar @nodes; $i++){

		
		my $sub1 = ${$subs_on_node{${$nodes[$i]}->get_name()}}{$site_index};
	#	my $derived1 = $myCodonTable -> translate($sub1->{"Substitution::derived_allele"});
		my $derived1 = ($sub1->{"Substitution::derived_allele"});
		my $ancestor1 = $sub1->{"Substitution::ancestral_allele"};
		
		my $different_usefulness = 0;
		my $prev_usefulness; 
			## stopped here
		my $count_same = 1; # to add the node1 itself
		my $count_diff;
print "node 1:  $ancestor1 $derived1 \n";
		for (my $j = 0; $j < scalar @nodes; $j++){
			if ($j == $i){ next; }
			if (value_is_in_array(${$nodes[$j]}, \@{ ${$nodes[$i]}->get_sisters })){ 
				print " Hello sister!\n";
				next; 
			}
			my $sub2 = ${$subs_on_node{${$nodes[$j]}->get_name()}}{$site_index};
	#		my $derived2 = $myCodonTable -> translate($sub2->{"Substitution::derived_allele"});
			my $derived2 = ($sub2->{"Substitution::derived_allele"});
			my $ancestor2 = $sub2->{"Substitution::ancestral_allele"};
			my $usefulness = $codon_evolution{$derived2} - $codon_evolution{$ancestor2};
			print " usefulness  $usefulness \n";
			if (defined $prev_usefulness){
				if ($prev_usefulness != $usefulness) {
					$different_usefulness = 1;
					print "Usefulness changed!\n";
				}
			}
			$prev_usefulness = $usefulness;
			if ($ancestor1 ne $ancestor2 ){ next; }
			print "node 2:  $ancestor2 $derived2 \n";
			my $dist = calc_my_distance(${$nodes[$i]}, ${$nodes[$j]});
			if ($dist > 0){ # ie these nodes are not sequential
				if (!exists $hash{"$ancestor1$derived1$i"} ){
					my @same = ();
					my @diff = ();
					$hash{"$ancestor1$derived1$i"} = (\@same, \@diff);
				}
				if ($derived1 eq $derived2){
					push @{ ($hash{"$ancestor1$derived1$i"})->[0] }, $dist;
					$count_same++;
				}
				else {
					push @{ ($hash{"$ancestor1$derived1$i"})->[1] }, $dist;
					$count_diff++;
				}
			}
			else {
				print "WOA! sequential nodes here!\n";
			}

		}
		
		print " count same: $count_same, count diff: $count_diff \n";
		push @{ ($hash{"$ancestor1$derived1$i"})->[2] }, $count_same;
		push @{ ($hash{"$ancestor1$derived1$i"})->[2] }, $count_diff;
		push @{ ($hash{"$ancestor1$derived1$i"})->[2] }, $different_usefulness; # for splitting into groups 

	}
	return %hash;		
}

#accepts sequentials (after 02 06 2015), does not accept sisters

sub shuffler {
	my $tree = $_[0];
	my @nodes = @{$_[1]};
	my %subs_on_node = %{$_[2]};
	my $site_index = $_[3];

	my %hash;
	my @distances_same;
	my @distances_diff;
	my $myCodonTable   = Bio::Tools::CodonTable->new();
	my %hash_of_nodes;
	
	#split the array of nodes by the ancestor codon 
	foreach my $node (@nodes){
		my $ancestor = ${$subs_on_node{${$node}->get_name()}}{$site_index}->{"Substitution::ancestral_allele"};
		push @{$hash_of_nodes{$ancestor}}, $node;
	}
	
	foreach my $ancestor (keys %hash_of_nodes){
		my @nodes_subset = @{$hash_of_nodes{$ancestor}};
		my @shuffled = shuffle @nodes_subset;

		for (my $i = 0; $i < scalar @nodes_subset; $i++){
			my $sub1 = ${$subs_on_node{${$nodes_subset[$i]}->get_name()}}{$site_index};
			my $derived1;
			if ($_[5] eq "nsyn"){
				$derived1 = $myCodonTable -> translate($sub1->{"Substitution::derived_allele"});
			}
			elsif ($_[5] eq "syn"){
				$derived1 = ($sub1->{"Substitution::derived_allele"});
			}
			else {
				die "wrong argument in sub shuffler; only syn or nsyn accepted as the 6th arg";
			}
			my $count_same = 1; # to add the node1 itself
			my $count_diff;
#print "node 1:  $ancestor $derived1 \n";
			for (my $j = 0; $j < scalar @nodes_subset; $j++){
				if ($j == $i){ next; }
				if (value_is_in_array(${$shuffled[$j]}, \@{ ${$shuffled[$i]}->get_sisters })){ #!
					next; 
				}
				my $sub2 = ${$subs_on_node{${$nodes_subset[$j]}->get_name()}}{$site_index}; ##mistake found: nodes instead of nodes_subset
				my $derived2;
				if ($_[5] eq "nsyn"){
					$derived2 = $myCodonTable -> translate($sub2->{"Substitution::derived_allele"});
				}
				elsif ($_[5] eq "syn"){
					$derived2 = ($sub2->{"Substitution::derived_allele"});
				}
				else {
					die "wrong argument in sub shuffler; only syn or nsyn accepted as the 6th arg";
				}
#print "node 2:   $ancestor $derived2 \n";
			#	my $dist = calc_my_distance(${$shuffled[$i]}, ${$shuffled[$j]});
			    my $dist = node_distance(${$shuffled[$i]}, ${$shuffled[$j]}); #!
#print " dist $dist\n";
				if ($dist > 0){ # ie these nodes are not sequential; does not work since 02 06 2015
					if (!exists $hash{"$ancestor$derived1$i"} ){
						my @same = ();
						my @diff = ();
						$hash{"$ancestor$derived1$i"} = (\@same, \@diff);
					}
					if ($derived1 eq $derived2){
						push @{ ($hash{"$ancestor$derived1$i"})->[0] }, $dist;
						$count_same++;
					}
					else {
						push @{ ($hash{"$ancestor$derived1$i"})->[1] }, $dist;
						$count_diff++;
					}
				}

		}
		
#print " count same: $count_same, count diff: $count_diff \n";
		push @{ ($hash{"$ancestor$derived1$i"})->[2] }, $count_same;
		push @{ ($hash{"$ancestor$derived1$i"})->[2] }, $count_diff;

	}
	}

	return %hash;		
}


## careful! it does not ignore sequentials, if method is calc_true_patristic_distance instead og calc_my_distance
sub find_all_distances_except_seq_and_sis_radius {
	my $tree = $_[0];
	my @nodes = @{$_[1]};
	my %subs_on_node = %{$_[2]};
	my $site_index = $_[3];

	my %hash;
	my @distances_same;
	my @distances_diff;
	my $myCodonTable   = Bio::Tools::CodonTable->new();
	
	for (my $i = 0; $i < scalar @nodes; $i++){

		
		my $sub1 = ${$subs_on_node{${$nodes[$i]}->get_name()}}{$site_index};
			my $derived1;
			if ($_[5] eq "nsyn"){
				$derived1 = $myCodonTable -> translate($sub1->{"Substitution::derived_allele"});
			}
			elsif ($_[5] eq "syn"){
				$derived1 = ($sub1->{"Substitution::derived_allele"});
			}
			else {
				die "wrong argument in sub shuffler; only syn or nsyn accepted as the 6th arg";
			}
		my $ancestor1 = $sub1->{"Substitution::ancestral_allele"};
		
		my $count_same = 1; # to add the node1 itself
		my $count_diff;
#print "node 1:  $ancestor1 $derived1 \n";
		for (my $j = 0; $j < scalar @nodes; $j++){
			if ($j == $i){ next; }
			if (value_is_in_array(${$nodes[$j]}, \@{ ${$nodes[$i]}->get_sisters })){ 
			#	print " Hello sister!\n";
				next; 
			}
			my $sub2 = ${$subs_on_node{${$nodes[$j]}->get_name()}}{$site_index};
			my $derived2;
			if ($_[5] eq "nsyn"){
				$derived2 = $myCodonTable -> translate($sub2->{"Substitution::derived_allele"});
			}
			elsif ($_[5] eq "syn"){
				$derived2 = ($sub2->{"Substitution::derived_allele"});
			}
			else {
				die "wrong argument in sub shuffler; only syn or nsyn accepted as the 6th arg";
			}
			my $ancestor2 = $sub2->{"Substitution::ancestral_allele"};

			if ($ancestor1 ne $ancestor2 ){ next; }
			#print "node 2:  $ancestor2 $derived2 \n";
			my $dist = calc_true_patristic_distance(${$nodes[$i]}, ${$nodes[$j]});
			if ($dist > 0){ # ie these nodes are not sequential
				if (!exists $hash{"$ancestor1$derived1$i"} ){
					my @same = ();
					my @diff = ();
					$hash{"$ancestor1$derived1$i"} = (\@same, \@diff);
				}
				if ($derived1 eq $derived2){
					push @{ ($hash{"$ancestor1$derived1$i"})->[0] }, $dist;
					$count_same++;
				}
				else {
					push @{ ($hash{"$ancestor1$derived1$i"})->[1] }, $dist;
					$count_diff++;
				}
			}
			else {
			#	print "WOA! sequential nodes here!\n";
			}

		}
		
#print " count same: $count_same, count diff: $count_diff \n";
		push @{ ($hash{"$ancestor1$derived1$i"})->[2] }, $count_same;
		push @{ ($hash{"$ancestor1$derived1$i"})->[2] }, $count_diff;

	}

	return %hash;		
}



## accepts sequentials, does not accept sisters

sub collector {
	my $tree = $_[0];
	my @nodes = @{$_[1]};
	my %subs_on_node = %{$_[2]};
	my $site_index = $_[3];

	my %hash;
	my @distances_same;
	my @distances_diff;
	my $myCodonTable   = Bio::Tools::CodonTable->new();
	
	for (my $i = 0; $i < scalar @nodes; $i++){

		my $sub1 = ${$subs_on_node{${$nodes[$i]}->get_name()}}{$site_index};
		my $derived1 = $myCodonTable -> translate($sub1->{"Substitution::derived_allele"});
	#	my $derived1 = ($sub1->{"Substitution::derived_allele"});
		my $ancestor1 = $sub1->{"Substitution::ancestral_allele"};
		
		my $count_same = 1; # to add the node1 itself
		my $count_diff;
print "node 1:  $ancestor1 $derived1 \n";
		for (my $j = $i+1; $j < scalar @nodes; $j++){
			if (value_is_in_array(${$nodes[$j]}, \@{ ${$nodes[$i]}->get_sisters })){ 
				print " Hello sister!\n";
				next; 
			}
			my $sub2 = ${$subs_on_node{${$nodes[$j]}->get_name()}}{$site_index};
			my $derived2 = $myCodonTable -> translate($sub2->{"Substitution::derived_allele"});
	#		my $derived2 = ($sub2->{"Substitution::derived_allele"});
			my $ancestor2 = $sub2->{"Substitution::ancestral_allele"};

			if ($ancestor1 ne $ancestor2 ){ next; }
			print "node 2:  $ancestor2 $derived2 \n";
			my $dist = calc_true_patristic_distance(${$nodes[$i]}, ${$nodes[$j]});

				if (!exists $hash{"$ancestor1"} ){
					my @same = ();
					my @diff = ();
					$hash{"$ancestor1"} = (\@same, \@diff);
				}
				if ($derived1 eq $derived2){
					print "pushed $dist into same \n";
					push @{ ($hash{"$ancestor1"})->[0] }, $dist;
					$count_same++;
				}
				else {
					print "pushed $dist into diff \n";
					push @{ ($hash{"$ancestor1"})->[1] }, $dist;
					$count_diff++;
				}
		

		}
		
		print " count same: $count_same, count diff: $count_diff \n";
		push @{ ($hash{"$ancestor1"})->[2] }, $count_same;
		push @{ ($hash{"$ancestor1"})->[2] }, $count_diff;

	}
	return %hash;		
}


# ignore sequential nodes; all counts are kept separately for each node
sub find_all_distances_except_sequential_radius {
	my $tree = $_[0];
	my @nodes = @{$_[1]};
	my %subs_on_node = %{$_[2]};
	my $site_index = $_[3];

	my %hash;
	my @distances_same;
	my @distances_diff;
	my $myCodonTable   = Bio::Tools::CodonTable->new();
	
	for (my $i = 0; $i < scalar @nodes; $i++){

		
		my $sub1 = ${$subs_on_node{${$nodes[$i]}->get_name()}}{$site_index};
	#	my $derived1 = $myCodonTable -> translate($sub1->{"Substitution::derived_allele"});
		my $derived1 = ($sub1->{"Substitution::derived_allele"});
		my $ancestor1 = $sub1->{"Substitution::ancestral_allele"};
			
		my $count_same = 1; # to add the node1 itself
		my $count_diff;
print "node 1:  $ancestor1 $derived1 \n";
		for (my $j = 0; $j < scalar @nodes; $j++){
			if ($j == $i){ next; }
			my $sub2 = ${$subs_on_node{${$nodes[$j]}->get_name()}}{$site_index};
		#	my $derived2 = $myCodonTable -> translate($sub2->{"Substitution::derived_allele"});
			my $derived2 = ($sub2->{"Substitution::derived_allele"});
			my $ancestor2 = $sub2->{"Substitution::ancestral_allele"};
			if ($ancestor1 ne $ancestor2 ){ next; }
			print "node 2:  $ancestor2 $derived2 \n";
			my $dist = calc_my_distance(${$nodes[$i]}, ${$nodes[$j]});
			if ($dist > 0){ # ie these nodes are not sequential
				if (!exists $hash{"$ancestor1$derived1$i"} ){
					my @same = ();
					my @diff = ();
					$hash{"$ancestor1$derived1$i"} = (\@same, \@diff);
				}
				if ($derived1 eq $derived2){
					push @{ ($hash{"$ancestor1$derived1$i"})->[0] }, $dist;
					$count_same++;
				}
				else {
					push @{ ($hash{"$ancestor1$derived1$i"})->[1] }, $dist;
					$count_diff++;
				}
			}
			else {
				print "WOA! sequential nodes here!\n";
			}

		}
		
		print " count same: $count_same, count diff: $count_diff \n";
		push @{ ($hash{"$ancestor1$derived1$i"})->[2] }, $count_same;
		push @{ ($hash{"$ancestor1$derived1$i"})->[2] }, $count_diff;
		

	}
	return %hash;	
}

## undone and abandoned
sub find_all_distances_gold {
	my $tree = $_[0];
	my @nodes = @{$_[1]};
	my %subs_on_node = %{$_[2]};
	my $site_index = $_[3];

	my @distances_same;
	my @distances_diff;
	
	for (my $i = 0; $i < scalar @nodes; $i++){
		
		my $sub1 = ${$subs_on_node{${$nodes[$i]}->get_name()}}{$site_index};
		my $derived1 = $sub1->{"Substitution::derived_allele"};
		my $ancestor1 = $sub1->{"Substitution::ancestral_allele"};
		
			
		for (my $j = 0; $j < scalar @nodes; $j++){
			if ($j == $i){ next; }
			my $sub2 = ${$subs_on_node{${$nodes[$j]}->get_name()}}{$site_index};
			my $derived2 = $sub2->{"Substitution::derived_allele"};
			my $ancestor2 = $sub2->{"Substitution::ancestral_allele"};
			if ($ancestor1 ne $ancestor2 ){ next; }
			my $dist = calc_true_patristic_distance(${$nodes[$i]}, ${$nodes[$j]});
			set_distance($nodes[$i], $nodes[$j], $dist);
			if ($derived1 eq $derived2){
				push @distances_same, $dist;
				convergent_plus($nodes[$i], $nodes[$j]);
			}
			else {
				push @distances_diff, $dist;
				divergent_plus($nodes[$i], $nodes[$j]);
			}
		}
	}
	return  (\@distances_same, \@distances_diff);	
}

## collects all distances,  ignores oaks between maples. Same ancestor. returns hash. For synmutmap, logic_3

sub find_all_distances_3 {
	my $tree = $_[0];
	my @nodes = @{$_[1]};
	my %subs_on_node = %{$_[2]};
	my $site_index = $_[3];
	#my $same_derived = $_[4];

    my %hash;
	my @distances_same;
	my @distances_diff;
	
	for (my $i = 0; $i < scalar @nodes; $i++){

#		my $min_dist_same;
#		my $min_dist_diff;
		
		my $sub1 = ${$subs_on_node{${$nodes[$i]}->get_name()}}{$site_index};
		my $derived1 = $sub1->{"Substitution::derived_allele"};
		my $ancestor1 = $sub1->{"Substitution::ancestral_allele"};
			
		for (my $j = 0; $j < scalar @nodes; $j++){
			if ($j == $i){ next; }
			my $sub2 = ${$subs_on_node{${$nodes[$j]}->get_name()}}{$site_index};
			my $derived2 = $sub2->{"Substitution::derived_allele"};
			my $ancestor2 = $sub2->{"Substitution::ancestral_allele"};
			if ($ancestor1 ne $ancestor2 ){ next; }
			my $dist = calc_true_patristic_distance(${$nodes[$i]}, ${$nodes[$j]});
			if (!exists $hash{$ancestor1} ){
				my @same = ();
				my @diff = ();
				$hash{$ancestor1} = (\@same, \@diff);
			}

			if ($derived1 eq $derived2){
				push @{ ($hash{$ancestor1})->[0] }, $dist;
			}
			else {
				push @{ ($hash{$ancestor1})->[1] }, $dist;
			}


		}

	}
	return %hash;	
}



## method for logic_general_same_ancestor: only considers mutations with the same ancestor allele
## supposed to be used with codonmutmap: tries to translate the resulting codon into an aa
sub find_min_distances_same_ancestor{
	my $tree = $_[0];
	my @nodes = @{$_[1]};
	my %subs_on_node = %{$_[2]};
	my $site_index = $_[3];
	#my $same_derived = $_[4];
	my $myCodonTable   = Bio::Tools::CodonTable->new();
	
	
	my @min_distances_same;
	my @min_distances_diff;
	
	for (my $i = 0; $i < scalar @nodes; $i++){
		my $min_dist_same;
		my $min_dist_diff;
		my $bit_vect1 = ${$nodes[$i]}->get_generic("-mutations_on_path")->Clone();
		$bit_vect1->Move_Left(1);
		my $sub1 = ${$subs_on_node{${$nodes[$i]}->get_name()}}{$site_index};
		my $derived1 = $myCodonTable->translate($sub1->{"Substitution::derived_allele"});
		my $ancestor1 = $sub1->{"Substitution::ancestral_allele"};

		for (my $j = 0; $j < scalar @nodes; $j++){
			if ($j == $i ){ next; }
			my $sub2 = ${$subs_on_node{${$nodes[$j]}->get_name()}}{$site_index};
			my $derived2 = $myCodonTable->translate($sub2->{"Substitution::derived_allele"});
			my $ancestor2 = $sub2->{"Substitution::ancestral_allele"};
			if ($ancestor1 ne $ancestor2 ){ next; }
			my $mrca = ${$nodes[$i]}->get_mrca(${$nodes[$j]})->get_generic("-mutations_on_path")->Size();
			my $bit_vect2 = ${$nodes[$j]}->get_generic("-mutations_on_path")->Clone();
			my $bit_vect1_t = $bit_vect1->Clone();
			$bit_vect2->Move_Left(1);
			$bit_vect2->Move_Right($mrca+1);
			$bit_vect1_t->Move_Right($mrca+1);
			if ($bit_vect2->is_empty() & $bit_vect1_t->is_empty()){
				my $dist = calc_true_patristic_distance(${$nodes[$i]}, ${$nodes[$j]});
				
				if ($derived1 eq $derived2){
					if(!defined $min_dist_same || $dist < $min_dist_same){
						$min_dist_same = $dist;
					}
				}
				else {
					if(!defined $min_dist_diff || $dist < $min_dist_diff){
						$min_dist_diff = $dist;
					}
				}
			}
		}
		push @min_distances_same, $min_dist_same;
		push @min_distances_diff, $min_dist_diff;
	}
	return (\@min_distances_same, \@min_distances_diff);	
}



## for usage with synmutmap, logic_3
sub find_min_distances_3{
	my $tree = $_[0];
	my @nodes = @{$_[1]};
	my %subs_on_node = %{$_[2]};
	my $site_index = $_[3];
	
	my %hash;
	
	for (my $i = 0; $i < scalar @nodes; $i++){
		my $min_dist_same;
		my $min_dist_diff;
		my $bit_vect1 = ${$nodes[$i]}->get_generic("-mutations_on_path")->Clone();
		$bit_vect1->Move_Left(1);
		my $sub1 = ${$subs_on_node{${$nodes[$i]}->get_name()}}{$site_index};
		my $derived1 = $sub1->{"Substitution::derived_allele"};
		my $ancestor1 = $sub1->{"Substitution::ancestral_allele"};

		for (my $j = 0; $j < scalar @nodes; $j++){
			if ($j == $i ){ next; }
			my $sub2 = ${$subs_on_node{${$nodes[$j]}->get_name()}}{$site_index};
			my $derived2 = $sub2->{"Substitution::derived_allele"};
			my $ancestor2 = $sub2->{"Substitution::ancestral_allele"};
			if ($ancestor1 ne $ancestor2 ){ next; }
			my $mrca = ${$nodes[$i]}->get_mrca(${$nodes[$j]})->get_generic("-mutations_on_path")->Size();
			my $bit_vect2 = ${$nodes[$j]}->get_generic("-mutations_on_path")->Clone();
			my $bit_vect1_t = $bit_vect1->Clone();
			$bit_vect2->Move_Left(1);
			$bit_vect2->Move_Right($mrca+1);
			$bit_vect1_t->Move_Right($mrca+1);
			if ($bit_vect2->is_empty() & $bit_vect1_t->is_empty()){
				my $dist = calc_true_patristic_distance(${$nodes[$i]}, ${$nodes[$j]});
				
				if ($derived1 eq $derived2){
					if(!defined $min_dist_same || $dist < $min_dist_same){
						$min_dist_same = $dist;
					}
				}
				else {
					if(!defined $min_dist_diff || $dist < $min_dist_diff){
						$min_dist_diff = $dist;
					}
				}
			}
		}
		
		if (!exists $hash{$ancestor1} ){
			my @same;
			my @diff;
			$hash{$ancestor1} = (\@same, \@diff);
		}

		push @{ ($hash{$ancestor1})->[0] }, $min_dist_same;
		push @{ ($hash{$ancestor1})->[1] }, $min_dist_diff;
		print scalar @{ ($hash{$ancestor1})->[0] }." pushing $min_dist_same, $min_dist_diff";
		print "\n";
	}
	return %hash;	
}





## for usage with synmutmap
sub find_min_distances_same_ancestor_syn{
	my $tree = $_[0];
	my @nodes = @{$_[1]};
	my %subs_on_node = %{$_[2]};
	my $site_index = $_[3];
	
	
	my @min_distances_same;
	my @min_distances_diff;
	
	for (my $i = 0; $i < scalar @nodes; $i++){
		my $min_dist_same;
		my $min_dist_diff;
		my $bit_vect1 = ${$nodes[$i]}->get_generic("-mutations_on_path")->Clone();
		$bit_vect1->Move_Left(1);
		my $sub1 = ${$subs_on_node{${$nodes[$i]}->get_name()}}{$site_index};
		my $derived1 = $sub1->{"Substitution::derived_allele"};
		my $ancestor1 = $sub1->{"Substitution::ancestral_allele"};

		for (my $j = 0; $j < scalar @nodes; $j++){
			if ($j == $i ){ next; }
			my $sub2 = ${$subs_on_node{${$nodes[$j]}->get_name()}}{$site_index};
			my $derived2 = $sub2->{"Substitution::derived_allele"};
			my $ancestor2 = $sub2->{"Substitution::ancestral_allele"};
			if ($ancestor1 ne $ancestor2 ){ next; }
			my $mrca = ${$nodes[$i]}->get_mrca(${$nodes[$j]})->get_generic("-mutations_on_path")->Size();
			my $bit_vect2 = ${$nodes[$j]}->get_generic("-mutations_on_path")->Clone();
			my $bit_vect1_t = $bit_vect1->Clone();
			$bit_vect2->Move_Left(1);
			$bit_vect2->Move_Right($mrca+1);
			$bit_vect1_t->Move_Right($mrca+1);
			if ($bit_vect2->is_empty() & $bit_vect1_t->is_empty()){
				my $dist = calc_true_patristic_distance(${$nodes[$i]}, ${$nodes[$j]});
				
				if ($derived1 eq $derived2){
					if(!defined $min_dist_same || $dist < $min_dist_same){
						$min_dist_same = $dist;
					}
				}
				else {
					if(!defined $min_dist_diff || $dist < $min_dist_diff){
						$min_dist_diff = $dist;
					}
				}
			}
		}
		push @min_distances_same, $min_dist_same;
		push @min_distances_diff, $min_dist_diff;
	}
	return (\@min_distances_same, \@min_distances_diff);	
}

sub find_min_distances_naive{
	my $tree = $_[0];
	my @nodes = @{$_[1]};
	my %subs_on_node = %{$_[2]};
	my $site_index = $_[3];
	#my $same_derived = $_[4];
	my @min_distances_same;
	my @min_distances_diff;
	
	for (my $i = 0; $i < scalar @nodes; $i++){
		my $min_dist_same;
		my $min_dist_diff;
		my $bit_vect1 = ${$nodes[$i]}->get_generic("-mutations_on_path")->Clone();
		$bit_vect1->Move_Left(1);
		my $derived1 = ${$subs_on_node{${$nodes[$i]}->get_name()}}{$site_index}->{"Substitution::derived_allele"};
		#my @sieved = sieve(\@nodes,\%subs_on_node, $site_index, $t_derived, $same_derived);
		#for (my $j = 0; $j < scalar @sieved; $j++){
		for (my $j = 0; $j < scalar @nodes; $j++){
			if ($j == $i){ next; }
			#my $mrca = ${$nodes[$i]}->get_mrca(${$sieved[$j]})->get_generic("-mutations_on_path")->Size();
			my $mrca = ${$nodes[$i]}->get_mrca(${$nodes[$j]})->get_generic("-mutations_on_path")->Size();
			#my $bit_vect2 = ${$sieved[$j]}->get_generic("-mutations_on_path")->Clone();
			my $bit_vect2 = ${$nodes[$j]}->get_generic("-mutations_on_path")->Clone();
			my $bit_vect1_t = $bit_vect1->Clone();
			$bit_vect2->Move_Left(1);
			$bit_vect2->Move_Right($mrca+1);
			$bit_vect1_t->Move_Right($mrca+1);
			if ($bit_vect2->is_empty() & $bit_vect1_t->is_empty()){
				#my $dist = ${$nodes[$i]}->calc_patristic_distance(${$sieved[$j]});
				my $dist = calc_true_patristic_distance(${$nodes[$i]}, ${$nodes[$j]});
				
				my $derived2 = ${$subs_on_node{${$nodes[$j]}->get_name()}}{$site_index}->{"Substitution::derived_allele"};
				if ($derived1 eq $derived2){
					if(!defined $min_dist_same || $dist < $min_dist_same){
						$min_dist_same = $dist;
					}
				}
				else {
					if(!defined $min_dist_diff || $dist < $min_dist_diff){
						$min_dist_diff = $dist;
					}
				}
			}
		}
		push @min_distances_same, $min_dist_same;
		push @min_distances_diff, $min_dist_diff;
	}
	return (\@min_distances_same, \@min_distances_diff);	
}



#test_my_visit_depth_first();

sub test_tree_lengths {
	my $tree = parse_tree("C:/Users/weidewind/Documents/CMD/Coevolution/Influenza/Kryazhimsky11/Mock/h1.l.r.newick");
	my $protein_name = "h1";
	tree_lengths($protein_name,$tree);
}

sub test_subtree_lengths {
	my $tree = parse_tree("C:/Users/weidewind/Documents/CMD/Coevolution/Influenza/Kryazhimsky11/Mock/h1.l.r.newick");
	my $protein_name = "h1";
	my $hash = subtree_lengths($protein_name, $tree);
	print $hash->{"alastar7"}." must be 13";
}

sub test_my_visit_depth_first {
	    set_mutmap("h1", "nsyn"); # for this sub to work properly path in set_mutmap had to be changed to Mock
		my $tree = parse_tree("C:/Users/weidewind/Documents/CMD/Coevolution/Influenza/Kryazhimsky11/Mock/h1.l.r.newick");
	    my @array;
	    my $root = $tree-> get_root;
	    my @args = (2,5, $root);
	    my_visit_depth_first ($root, \@array,\&update_ring,\&has_no_mutation,\@args,0);
	    print "1: up to 5 ".$static_ring_hash{2}{"alister1"}{1}[0]."\t_".$static_ring_hash{2}{"alister1"}{1}[1]."\n";
	    print "2: up to 10 ".$static_ring_hash{2}{"alister1"}{2}[0]."\t_".$static_ring_hash{2}{"alister1"}{2}[1]."\n";
	    
}


my @h1uptrend = qw(111 169 205 113);
#entrenchment_bootstrap("h1", 100, 150);
#boot_median_test("h1", 150);

#set_mutmap("h3", "nsyn");
#set_distance_matrix("h3");
#depth_groups_entrenchment_optimized(10,150, \@h3_internal);

# "green" control
#entrenchment_bootstrap_full("n2", 100, 150);




#entrenchment_blue_violet_bootstrap("n1", 100, 0); #restriction - last


#enrichment_blue_violet_optimized ("C:/Users/weidewind/workspace/perlCoevolution/TreeUtils/Phylo/MutMap/n2_for_enrichment_1_norestr_obs", "C:/Users/weidewind/workspace/perlCoevolution/TreeUtils/Phylo/MutMap/n2_for_enrichment_1_norestr", \@n2_internal);

#enrichment_optimized("C:/Users/weidewind/workspace/perlCoevolution/TreeUtils/Phylo/MutMap/n1_opti_proc_100_group_obs", "C:/Users/weidewind/workspace/perlCoevolution/TreeUtils/Phylo/MutMap/n1_opti_proc_100_group", \@n1_wan_epitopes);

#write_csv("h1");
#write_csv("h3");
#write_csv("n1");
#write_csv("n2");

#rewrite_csv("h3");
#set_mutmap("h1", "nsyn");
#boot_median_test("h1");

sub entrenchment_bootstrap{ 
	my $prot = $_[0];
	my $restriction = $_[1];
	set_mutmap($prot, "nsyn");
	set_distance_matrix($prot);
	my %matrix = incidence_matrix(); #!
	print_incidence_matrix(\%matrix, "C:/Users/weidewind/workspace/perlCoevolution/TreeUtils/Phylo/MutMap/");
	my %hash = depth_groups_entrenchment_optimized(10,$restriction);
	my $norm;
	foreach my $bin(1..32){
			$norm += $hash{$bin}[0];
	}

	my @simulated_hists;
	for (my $i = 1; $i < 101; $i++){
		my @mock_mutmaps = read_incidence_matrix("C:/Users/weidewind/workspace/perlCoevolution/TreeUtils/Phylo/MutMap/".$prot."_ShuffledMatrices/$i");
		
		%static_subs_on_node = %{$mock_mutmaps[0]};
		%static_nodes_with_sub = %{$mock_mutmaps[1]};

		my %hash = depth_groups_entrenchment_optimized(10,$restriction);
		push @simulated_hists, \%hash;
		%static_ring_hash = ();
		%static_depth_hash = ();
	}

	store \@simulated_hists, "C:/Users/weidewind/workspace/perlCoevolution/TreeUtils/Phylo/MutMap/".$prot."_full_simulation_".$restriction."_stored";
	my $arref = retrieve("C:/Users/weidewind/workspace/perlCoevolution/TreeUtils/Phylo/MutMap/".$prot."_full_simulation_".$restriction."_stored");
	
	open CSV, ">C:/Users/weidewind/workspace/perlCoevolution/TreeUtils/Phylo/MutMap/".$prot."_full_simulation_".$restriction.".csv";
	foreach my $bin(1..32){
			print CSV $bin."_obs,".$bin."_exp,";
	}
	print CSV "\n";

	foreach my $i(0..99){
		my $sum;
		foreach my $bin(1..32){
				$sum += $arref->[$i]->{$bin}->[0];
		}
		foreach my $bin(1..32){
			print CSV ($arref->[$i]->{$bin}->[0])*$norm/$sum.",".($arref->[$i]->{$bin}->[1])*$norm/$sum.",";
		}
		print CSV"\n";
	}
	close CSV;
}

sub entrenchment_bootstrap_full{ 
	my $prot = $_[0];
	my $iterations = $_[1];
	my $restriction = $_[2];
	set_mutmap($prot, "nsyn");
	set_distance_matrix($prot);
	my %matrix = incidence_matrix(); #!
	print_incidence_matrix(\%matrix, "C:/Users/weidewind/workspace/perlCoevolution/TreeUtils/Phylo/MutMap/");
	my %obs_hash = depth_groups_entrenchment_optimized(10,$restriction);
	my $norm;
	foreach my $bin(1..32){
			$norm += $obs_hash{$bin}[0];
	}
	
	%static_ring_hash = ();
	%static_depth_hash = ();
	%static_subtree_info = ();
	
	my @simulated_hists;
	for (my $i = 1; $i <= $iterations; $i++){
		my @mock_mutmaps = read_incidence_matrix("C:/Users/weidewind/workspace/perlCoevolution/TreeUtils/Phylo/MutMap/".$prot."_ShuffledMatrices/$i");
		
		%static_subs_on_node = %{$mock_mutmaps[0]};
		%static_nodes_with_sub = %{$mock_mutmaps[1]};

		my %hash = depth_groups_entrenchment_optimized(10,$restriction);
		push @simulated_hists, \%hash;
		
		my $sum;
		foreach my $bin(1..32){
				$sum += $hash{$bin}[0];
		}
		foreach my $bin(1..32){
			$hash{$bin}[0] = $hash{$bin}[0]*$norm/$sum;
			$hash{$bin}[1] = $hash{$bin}[1]*$norm/$sum;
		}
		
		%static_ring_hash = ();
		%static_depth_hash = ();
		%static_subtree_info = ();
	}

	store \@simulated_hists, "C:/Users/weidewind/workspace/perlCoevolution/TreeUtils/Phylo/MutMap/".$prot."_noneed_simulation_".$restriction."_stored";
	my $arref = retrieve("C:/Users/weidewind/workspace/perlCoevolution/TreeUtils/Phylo/MutMap/".$prot."_noneed_simulation_".$restriction."_stored");
	
	open CSV, ">C:/Users/weidewind/workspace/perlCoevolution/TreeUtils/Phylo/MutMap/".$prot."_noneed_simulation_".$restriction.".csv";
	foreach my $bin(1..32){
			print CSV $bin."_obs,".$bin."_exp,";
	}
	print CSV "\n";

	foreach my $i(0..$iterations-1){
		foreach my $bin(1..32){
			print CSV ($arref->[$i]->{$bin}->[0]).",".($arref->[$i]->{$bin}->[1]).",";
		}
		print CSV"\n";
	}
	close CSV;
	
	my $file = $prot."_noneed_boot_median_test_".$restriction;
	open FILE, ">$file";
	foreach my $bin(1..32){
			$obs_hash{$bin} = $obs_hash{$bin}->[0];
		}
	my $obs_median = hist_median_for_hash(\%obs_hash);
	
	print FILE "\n observed median: $obs_median\n";
	
	my $pval_epi;
	my $pval_env;
	
	foreach my $i(0..$iterations-1){
		my %hash;
		foreach my $bin(1..32){
			$hash{$bin} = $arref->[$i]->{$bin}->[0];
		}
		my $boot_median = hist_median_for_hash(\%hash);
		print FILE "\n boot median: $boot_median\n";
		if ($boot_median <= $obs_median){
			$pval_epi += 1;
		}
		if ($boot_median >= $obs_median){
			$pval_env += 1;
		}

	}
	print FILE "pvalue epistasis ".($pval_epi/$iterations)." pvalue environment ".($pval_env/$iterations);
	close FILE;
}


#check_entrenchment_blue_violet_bootstrap ("h3", 200);
#	check_part_2("h3", 20, 150);

sub test_remove_zero_branches {
	my $tree = parse_tree("C:/Users/weidewind/Documents/CMD/Coevolution/Influenza/Kryazhimsky11/h1.l.r.newick");
	my $multitree = remove_zero_branches($tree);
	open TREE, ">multitree.tre" or die "Cannot create file";
	my $treestr=tree2str($multitree);
	print TREE $treestr;
	close TREE;
	
}

test_remove_zero_branches();

sub check_entrenchment_blue_violet_bootstrap {
	my $prot = $_[0];
	my $iterations = $_[1];
	set_mutmap($prot, "nsyn");
	set_distance_matrix($prot);
	my %matrix = incidence_matrix(); #!
	print_incidence_matrix(\%matrix, "C:/Users/weidewind/workspace/perlCoevolution/TreeUtils/Phylo/MutMap/");
	my %obs_hash = depth_groups_entrenchment_optimized(1,0);
	my %real_obs_hash;
	my %real_exp_hash;
	my $norm;
	
	my $maxbin = 0;
	
	foreach my $bin(keys %obs_hash){
			$maxbin = max($bin, $maxbin);
			$norm += $obs_hash{$bin}[0];
	}
	
	%static_ring_hash = ();
	%static_depth_hash = ();
	%static_subtree_info = ();
	
	my @simulated_hists;
	my %simulated_medians;
	my %all_simulated_medians;
	my %all_hash;
	
	open STORAGE, ">C:/Users/weidewind/workspace/perlCoevolution/TreeUtils/Phylo/MutMap/".$prot."_simply_norestr_stored_200.csv";
	for (my $i = 1; $i <= $iterations; $i++){
		my @mock_mutmaps = read_incidence_matrix("C:/Users/weidewind/workspace/perlCoevolution/TreeUtils/Phylo/MutMap/".$prot."_ShuffledMatrices/$i");
		
		%static_subs_on_node = %{$mock_mutmaps[0]};
		%static_nodes_with_sub = %{$mock_mutmaps[1]};

		my %hash = check_depth_groups_entrenchment_optimized(1);
		push @simulated_hists, \%hash;
		
		my $sum;
		
		foreach my $bin(1..$maxbin){
				foreach my $mutnum (keys %hash){
					if ($mutnum > 0){
						foreach my $site_node(keys %{$hash{$mutnum}}){
							$sum += $hash{$mutnum}{$bin}[0];
						}
					}	
			}
		}
		
				foreach my $mutnum (keys %hash){
					if ($mutnum > 0){
					#print STORAGE ">$mutnum,".$hash{$mutnum}{"counter"}."\n";
					foreach my $bin(1..$maxbin){
						$hash{$mutnum}{$bin}[0] = $hash{$mutnum}{$bin}[0]*$norm/$sum;
						$all_hash{$mutnum}{$bin}[0] += $hash{$mutnum}{$bin}[0];
						$hash{$mutnum}{$bin}[1] = $hash{$mutnum}{$bin}[1]*$norm/$sum;
						$all_hash{$mutnum}{$bin}[1] += $hash{$mutnum}{$bin}[1];
					#	print STORAGE "$bin,".$hash{$mutnum}{$bin}[0].",".$hash{$mutnum}{$bin}[1]."\n";
					}
					$all_hash{$mutnum}{"counter"} += $hash{$mutnum}{"counter"};
					my $boot_obs_median = hist_median_for_hash_arr(\%{$hash{$mutnum}}, 0);
					push @{$simulated_medians{$mutnum}[0]}, $boot_obs_median;
					my $boot_exp_median = hist_median_for_hash_arr(\%{$hash{$mutnum}}, 1);
					push @{$simulated_medians{$mutnum}[1]}, $boot_exp_median;
					}
			}
		
		%static_ring_hash = ();
		%static_depth_hash = ();
		%static_subtree_info = ();
	}
	
	my @sorted_keys = sort {$a <=> $b} keys %simulated_medians;
	my %total_medians;
	foreach my $mutnum(@sorted_keys){
		my $stat_obs = Statistics::Descriptive::Full->new();
		$stat_obs->add_data(\@{$simulated_medians{$mutnum}[0]});
		my $stat_exp = Statistics::Descriptive::Full->new();
		$stat_exp->add_data(\@{$simulated_medians{$mutnum}[1]});
		print STORAGE $mutnum.",".hist_median_for_hash_arr(\%{$all_hash{$mutnum}}, 0).",".
		hist_mean_for_hash_arr(\%{$all_hash{$mutnum}}, 0).",".stddev(\@{$simulated_medians{$mutnum}[0]}).",".
		$stat_obs->percentile(5).",".
		$stat_obs->percentile(95).",".
		hist_median_for_hash_arr(\%{$all_hash{$mutnum}}, 1).",".
		hist_mean_for_hash_arr(\%{$all_hash{$mutnum}}, 1).",".stddev(\@{$simulated_medians{$mutnum}[1]}).",".
		$stat_exp->percentile(5).",".
		$stat_exp->percentile(95).",".
		$all_hash{$mutnum}{"counter"}."\n";
	}
	
	
	
	foreach my $mutnum(@sorted_keys){
		print STORAGE $mutnum."\n";
		for (my $i = 1; $i <= $iterations; $i++){
			print STORAGE $simulated_medians{$mutnum}[0][$i].",".$simulated_medians{$mutnum}[1][$i]."\n";
		}
	}
close STORAGE;
	#store \@simulated_hists, "C:/Users/weidewind/workspace/perlCoevolution/TreeUtils/Phylo/MutMap/".$prot."_check_simulation_150_hash_stored";
	
}

sub check_part_2	{
		my $prot = $_[0];
	my $iterations = $_[1];
	my $restriction = $_[2];

	my $maxbin = 400;
		my $arref = retrieve("C:/Users/weidewind/workspace/perlCoevolution/TreeUtils/Phylo/MutMap/".$prot."_check_simulation_".$restriction."_stored");
	
	my $file = $prot."_check_statistics_".$restriction;
	open FILE, ">$file";
	
	my $pval_epi;
	my $pval_env;
	foreach my $mutnum (1..50){
		print $mutnum."\n";
			my %complete_boot_obs_hash;
			my %complete_boot_exp_hash;
			my $complete_counter;
		foreach my $i(0..$iterations-1){
			my %boot_obs_hash;
			my %boot_exp_hash;
			
			my $counter;
			foreach my $maxdepth (0..7){
				#print FILE "\t maxdepth $maxdepth\n";
				if ($maxdepth >= $restriction/50){
					#print FILE "maxdepth ok\n";
					if ($arref->[$i]->[$maxdepth]){
						foreach my $testkey (keys %{$arref->[$i]->[$maxdepth]}){
							#print " found mutnum $testkey\n";
						}
						foreach my $site_node (keys %{$arref->[$i]->[$maxdepth]->{$mutnum}}){
							#print FILE " \t site_node $site_node\n";
								$counter++;
								$complete_counter++;
								foreach my $bin(1..$maxbin){
									$boot_obs_hash{$bin} += $arref->[$i]->[$maxdepth]->{$mutnum}->{$site_node}->{$bin}->[0];
									$boot_exp_hash{$bin} += $arref->[$i]->[$maxdepth]->{$mutnum}->{$site_node}->{$bin}->[1];
									$complete_boot_obs_hash{$bin} += $arref->[$i]->[$maxdepth]->{$mutnum}->{$site_node}->{$bin}->[0];
									$complete_boot_exp_hash{$bin} += $arref->[$i]->[$maxdepth]->{$mutnum}->{$site_node}->{$bin}->[1];
									#print FILE " added".$arref->[$i]->[$maxdepth]->{$mutnum}->{$site_node}->{$bin}->[0]." and ".$arref->[$i]->{$maxdepth}->{$mutnum}->{$site_node}->{$bin}->[1]."\n";
								
								}
						
					}
					}
				}
			}
			
			
			my $boot_obs_median = hist_median_for_hash(\%boot_obs_hash);
			my $boot_exp_median = hist_median_for_hash(\%boot_exp_hash);
			
			#print FILE "\n boot medians: $boot_obs_median , $boot_exp_median".."\n";

		}
		my $complete_boot_obs_median = hist_median_for_hash(\%complete_boot_obs_hash);
		my $complete_boot_exp_median = hist_median_for_hash(\%complete_boot_exp_hash);
		print FILE "mutnum $mutnum obs median $complete_boot_obs_median, exp median $complete_boot_exp_median, counter $complete_counter\n";
	}

	close FILE;
}


sub entrenchment_blue_violet_bootstrap{ 
	my $prot = $_[0];
	my $iterations = $_[1];
	my $restriction = $_[2];
	set_mutmap($prot, "nsyn");
	set_distance_matrix($prot);
	my %matrix = incidence_matrix(); #!
	print_incidence_matrix(\%matrix, "C:/Users/weidewind/workspace/perlCoevolution/TreeUtils/Phylo/MutMap/");
	my %obs_hash = depth_groups_entrenchment_optimized(1,$restriction);
	my %real_obs_hash;
	my %real_exp_hash;
	my $norm;
	
	my $maxbin = 0;
	
	foreach my $bin(keys %obs_hash){
			$maxbin = max($bin, $maxbin);
			$norm += $obs_hash{$bin}[0];
	}
	
	%static_ring_hash = ();
	%static_depth_hash = ();
	%static_subtree_info = ();
	
	my @simulated_hists;
	for (my $i = 1; $i <= $iterations; $i++){
		my @mock_mutmaps = read_incidence_matrix("C:/Users/weidewind/workspace/perlCoevolution/TreeUtils/Phylo/MutMap/".$prot."_ShuffledMatrices/$i");
		
		%static_subs_on_node = %{$mock_mutmaps[0]};
		%static_nodes_with_sub = %{$mock_mutmaps[1]};

		my %hash = depth_groups_entrenchment_optimized(1,$restriction);
		push @simulated_hists, \%hash;
		
		my $sum;
		foreach my $bin(1..$maxbin){
				$sum += $hash{$bin}[0];
		}
		foreach my $bin(1..$maxbin){
			$hash{$bin}[0] = $hash{$bin}[0]*$norm/$sum;
			$hash{$bin}[1] = $hash{$bin}[1]*$norm/$sum;
		}
		
		%static_ring_hash = ();
		%static_depth_hash = ();
		%static_subtree_info = ();
	}

	store \@simulated_hists, "C:/Users/weidewind/workspace/perlCoevolution/TreeUtils/Phylo/MutMap/".$prot."_1_100noneed_bv_simulation_".$restriction."_stored";
	my $arref = retrieve("C:/Users/weidewind/workspace/perlCoevolution/TreeUtils/Phylo/MutMap/".$prot."_1_100noneed_bv_simulation_".$restriction."_stored");
	
	open CSV, ">C:/Users/weidewind/workspace/perlCoevolution/TreeUtils/Phylo/MutMap/".$prot."_1_100noneed_bv_simulation_".$restriction.".csv";
	foreach my $bin(1..$maxbin){
			print CSV $bin."_obs,".$bin."_exp,";
	}
	print CSV "\n";

	foreach my $i(0..$iterations-1){
		foreach my $bin(1..$maxbin){
			print CSV ($arref->[$i]->{$bin}->[0]).",".($arref->[$i]->{$bin}->[1]).",";
		}
		print CSV"\n";
	}
	close CSV;
	
	my $file = $prot."_1_100noneed_bv_boot_median_test_".$restriction;
	open FILE, ">$file";

	foreach my $bin(1..$maxbin){
			$real_obs_hash{$bin} = $obs_hash{$bin}->[0];
			$real_exp_hash{$bin} = $obs_hash{$bin}->[1];
		}
	my $real_obs_median = hist_median_for_hash(\%real_obs_hash);
	my $real_exp_median = hist_median_for_hash(\%real_exp_hash);
	
	print FILE "\n observed medians: $real_obs_median , $real_exp_median\n";
	
	my $pval_epi;
	my $pval_env;
	
	foreach my $i(0..$iterations-1){
		my %boot_obs_hash;
		my %boot_exp_hash;
		foreach my $bin(1..$maxbin){
			$boot_obs_hash{$bin} = $arref->[$i]->{$bin}->[0];
			$boot_exp_hash{$bin} = $arref->[$i]->{$bin}->[1];
		}
		my $boot_obs_median = hist_median_for_hash(\%boot_obs_hash);
		my $boot_exp_median = hist_median_for_hash(\%boot_exp_hash);
		print FILE "\n boot medians: $boot_obs_median , $boot_exp_median\n";
		if ($boot_obs_median - $boot_exp_median >= $real_obs_median - $real_exp_median){
			$pval_env += 1;
		}
		if ($boot_obs_median - $boot_exp_median <= $real_obs_median - $real_exp_median){
			$pval_epi += 1;
		}

	}
	print FILE "pvalue epistasis ".($pval_epi/$iterations)." pvalue environment ".($pval_env/$iterations);
	close FILE;
}


sub rewrite_csv {
	my $prot = $_[0];
	my $arref = retrieve("C:/Users/weidewind/workspace/perlCoevolution/TreeUtils/Phylo/MutMap/".$prot."_simulation_150_stored");
	open CSV, ">C:/Users/weidewind/workspace/perlCoevolution/TreeUtils/Phylo/MutMap/".$prot."_simulation_150.csv";
	
	foreach my $bin(1..32){
			print CSV $bin."_obs,".$bin."_exp,";
	}
	print CSV "\n";
	foreach my $i(0..99){
		my $sum;
		foreach my $bin(1..32){
				$sum += $arref->[$i]->{$bin}->[0];
		}
		foreach my $bin(1..32){
			print CSV ($arref->[$i]->{$bin}->[0])*326/$sum.",".($arref->[$i]->{$bin}->[1])*326/$sum.",";
		}
		print CSV"\n";
	}
	close CSV;
	
}


sub write_csv {
	my $prot = $_[0];
	my $arref = retrieve("C:/Users/weidewind/workspace/perlCoevolution/TreeUtils/Phylo/MutMap/".$prot."_simulation_50_stored");
	open CSV, ">C:/Users/weidewind/workspace/perlCoevolution/TreeUtils/Phylo/MutMap/".$prot."_simulation_50.csv";
	
	foreach my $bin(1..31){
			print CSV $bin."_obs,".$bin."_exp,";
	}
	print CSV "\n";

	foreach my $i(0..101){
		foreach my $bin(1..31){
			print CSV $arref->[$i]->{$bin}->[0].",".$arref->[$i]->{$bin}->[1].",";
		}
		print CSV "\n";
	}
	close CSV;
}
#print "tester ".$arref->[0]->{3}->[0]."\t";
#print "tester ".$arref->[1]."\t";

sub main_optim {
	set_mutmap("h1", "nsyn");
	set_distance_matrix("h1");
	depth_groups_entrenchment_optimized(10);
}




#depth_groups_entrenchment(10);
#egor_site_entrenchment();
#global_entrenchment_epsilon(3, 1);
#sites_num();
# neva_site_entrenchment(10);
#depth_groups_entrenchment(10, \@h1uptrend);
#depth_groups_entrenchment_bootstrap(1, \@h3_antigenic_koel, 100);
 #maxdepths_hist();
 
sub boot_median_test {
	my $prot = $_[0];
	my $restriction = $_[1];
	my $arref = retrieve("C:/Users/weidewind/workspace/perlCoevolution/TreeUtils/Phylo/MutMap/".$prot."_simulation_".$restriction."_stored");

	my %obs_hash = depth_groups_entrenchment(10, $restriction);
	my $file = $prot."_new_boot_median_test_".$restriction;
	open FILE, ">$file";
	foreach my $bin(1..32){
			$obs_hash{$bin} = $obs_hash{$bin}->[0];
		}
	
	my $obs_median = hist_median_for_hash(\%obs_hash);
	
	print FILE "\n observed median: $obs_median\n";
	
	my $pval_epi;
	my $pval_env;
	
	
	foreach my $i(0..99){
		my %hash;
		my $sum;
		foreach my $bin(1..32){
				$sum += $arref->[$i]->{$bin}->[0];
		}
		foreach my $bin(1..32){
			$hash{$bin} = ($arref->[$i]->{$bin}->[0])*326/$sum;
		}
		my $boot_median = hist_median_for_hash(\%hash);
		print FILE "\n boot median: $boot_median\n";
		if ($boot_median <= $obs_median){
			$pval_epi += 1;
		}
		if ($boot_median >= $obs_median){
			$pval_env += 1;
		}

	}
	
	print FILE "pvalue epistasis ".($pval_epi/100)." pvalue environment ".($pval_env/100);
	close FILE;
} 
 
sub maxdepths_hist {
	my @array;
	my %hist;

	
	foreach my $ind (1..565){
		foreach my $node(@{$static_nodes_with_sub{$ind}}){
			$node = ${$node};
			if (!$node->is_terminal){
			#print $node->get_name();
			my @args = ($ind, 1, $node);
			my_visit_depth_first ($node, \@array,\&max_depth,\&has_no_mutation,\@args,0);
			print $static_depth_hash{$ind}{$node->get_name()}."\t";
			my_visit_depth_first ($node, \@array,\&update_ring,\&has_no_mutation,\@args,0);
			my $total_muts;
			my $total_length;
			foreach my $bin (keys %{$static_ring_hash{$ind}{$node->get_name()}}){
				$total_muts += $static_ring_hash{$ind}{$node->get_name()}{$bin}[0];
				$total_length += $static_ring_hash{$ind}{$node->get_name()}{$bin}[1];
			}
			if ($total_length == 0){
				print "NA\n"
			}
			else {
				print $total_muts/$total_length."\n";
			}
		}
		}
	}
} 
 
 
sub global_entrenchment {
	my $step = $_[0];
	my @group;
	if ($_[1]){
		@group = @{$_[1]};
	}
	else {
		@group = (1..565);
	}
	my $root = $static_tree-> get_root;
	my @array;
	my %hist;

	print "bin,observed,expected\n";
	foreach my $ind (@group){
		foreach my $node(@{$static_nodes_with_sub{$ind}}){
			$node = ${$node};
			if (!$node->is_terminal){
			#print $node->get_name();
			my @args = ($ind, $step, $node);
			my_visit_depth_first ($node, \@array,\&update_ring,\&has_no_mutation,\@args,0);
			my $total_muts;
			my $total_length;
			foreach my $bin (keys %{$static_ring_hash{$ind}{$node->get_name()}}){
				$total_muts += $static_ring_hash{$ind}{$node->get_name()}{$bin}[0];
				$total_length += $static_ring_hash{$ind}{$node->get_name()}{$bin}[1];
			}
			foreach my $bin (keys %{$static_ring_hash{$ind}{$node->get_name()}}){
				#print "bin $bin observed ".$static_ring_hash{$ind}{$node->get_name()}{$bin}[0]." totmut $total_muts totlen $total_length\n";
				if ($total_length > 0){ #there are some internal nodes with 0-length terminal daughter branches
					$hist{$bin}[0] += $static_ring_hash{$ind}{$node->get_name()}{$bin}[0]; #observed
					$hist{$bin}[1] += $total_muts*$static_ring_hash{$ind}{$node->get_name()}{$bin}[1]/$total_length; #expected
				}
			}
		}
		}
	}
	
	foreach my $bin (keys %hist){
		print $bin*$step.",".$hist{$bin}[0].",".$hist{$bin}[1]."\n";
	}
	
	
}


## full ring; each subtree has its own epsilon
sub global_entrenchment_epsilon {
	my $bin_count = $_[0];
	my @group;
	if ($_[2]){
		@group = @{$_[2]};
	}
	else {
		@group = (1..565);
	}
	my $scaled = $_[1];
	my $root = $static_tree-> get_root;
	my @array;
	my %hist;

	print "bin,observed,expected\n";
	foreach my $ind (@group){
		foreach my $node(@{$static_nodes_with_sub{$ind}}){
			$node = ${$node};
			if (!$node->is_terminal){
			print $node->get_name();
			
			my @args = ($ind, 1, $node);
			my_visit_depth_first ($node, \@array,\&max_depth,\&has_no_mutation,\@args,0);
			my $step = $static_depth_hash{$ind}{$node->get_name()}/$bin_count;
			print " maxdepth ".$static_depth_hash{$ind}{$node->get_name()}." step ".$step."\n";
			if ($step == 0){
				next;
			}
			@args = ($ind, $step, $node);
			my_visit_depth_first ($node, \@array,\&update_ring,\&has_no_mutation,\@args,0);
			
			#if ($static_depth_hash{$ind}{$node->get_name()} < 100){
			my $total_muts;
			my $total_length;
			foreach my $bin (keys %{$static_ring_hash{$ind}{$node->get_name()}}){
				$total_muts += $static_ring_hash{$ind}{$node->get_name()}{$bin}[0];
				$total_length += $static_ring_hash{$ind}{$node->get_name()}{$bin}[1];
			}
			#my $cumulative_muts;
			#my $cumulative_length;
			foreach my $bin (sort {$a <=> $b} keys %{$static_ring_hash{$ind}{$node->get_name()}}){
				my $histbin = $bin;
				if (!$scaled){
					$histbin = $step*$bin; # real radius
				}
				
				print "bin $bin observed ".$static_ring_hash{$ind}{$node->get_name()}{$bin}[0]." totmut $total_muts totlen $total_length\n";
				if ($total_length > 0){ #there are some internal nodes with 0-length terminal daughter branches
					#$cumulative_muts += $static_ring_hash{$ind}{$node->get_name()}{$bin}[0];
					#$cumulative_length += $static_ring_hash{$ind}{$node->get_name()}{$bin}[1];
					#$hist{$histbin}[0] += $cumulative_muts; #observed
					#$hist{$histbin}[1] += $total_muts*$cumulative_length/$total_length; #expected
				$hist{$histbin}[0] +=  $static_ring_hash{$ind}{$node->get_name()}{$bin}[0]; #observed
				$hist{$histbin}[1] += $total_muts*$static_ring_hash{$ind}{$node->get_name()}{$bin}[1]/$total_length;
				}
			}
		}
		}
	}
	
	foreach my $bin (sort {$a <=> $b} keys %hist){
		print $bin.",".$hist{$bin}[0].",".$hist{$bin}[1]."\n";
	}
	
	
}


sub find_epsilon {
	my @group;
	if ($_[1]){
		@group = @{$_[1]};
	}
	else {
		@group = (1..566);
	}
	my $root = $static_tree-> get_root;
	my @array;
	my @lambdas;

	foreach my $ind (@group){
		foreach my $node(@{$static_nodes_with_sub{$ind}}){
			$node = ${$node};
			if (!$node->is_terminal){
			my @args = ($ind, 600, $node);
			my_visit_depth_first ($node, \@array,\&update_ring,\&has_no_mutation,\@args,0);
			my $total_muts = $static_ring_hash{$ind}{$node->get_name()}{1}[0];
			if ($total_muts > 0){
				my $total_length = $static_ring_hash{$ind}{$node->get_name()}{1}[1];
				my $lambda = $total_muts;
				push @lambdas, $lambda;
				print $lambda."\n";
			}
		}
		}
	}
	
	return median(\@lambdas);

}


sub neva_site_entrenchment {
	my $step = $_[0];
	my $root = $static_tree-> get_root;
	my @array;
	my %hist;
	print "radius,site,node,observed,expected\n";
	my @group;
	if ($_[1]){
		@group = @{$_[1]};
	}
	else {
		@group = (1..565);
	}
	foreach my $ind (@group){
		foreach my $node(@{$static_nodes_with_sub{$ind}}){
			$node = ${$node};
			if (!$node->is_terminal){
			my @args = ($ind, $step, $node);
			my_visit_depth_first ($node, \@array,\&update_ring,\&has_no_mutation,\@args,0);
			my $total_muts;
			my $total_length;
			foreach my $bin (keys %{$static_ring_hash{$ind}{$node->get_name()}}){
				$total_muts += $static_ring_hash{$ind}{$node->get_name()}{$bin}[0];
				$total_length += $static_ring_hash{$ind}{$node->get_name()}{$bin}[1];
			}
			foreach my $bin (sort {$a <=> $b} (keys %{$static_ring_hash{$ind}{$node->get_name()}})){
				if ($total_length > 0 && $static_ring_hash{$ind}{$node->get_name()}{$bin}[1] > 0){ #there are some internal nodes with 0-length terminal daughter branches
					$hist{$bin}[0] = $static_ring_hash{$ind}{$node->get_name()}{$bin}[0]/$static_ring_hash{$ind}{$node->get_name()}{$bin}[1]; #observed
					$hist{$bin}[1] = $total_muts/$total_length; #expected
				}
				if (!$hist{$bin}[0]){
					$hist{$bin}[0] = 0;
				}
				if (!$hist{$bin}[1]){
					$hist{$bin}[1] = 0;
				}
				print "$bin,$ind,".$node->get_name().",".$hist{$bin}[0].",".$hist{$bin}[1]."\n";
			}
		}
		}
	}
	
	foreach my $bin (keys %hist){
		print " up to ".$bin*$step."\t".$hist{$bin}[0]."\t".$hist{$bin}[1]."\n";
	}
	
	
}

#bullshit_test();
sub bullshit_test{
	
	set_mutmap("h1", "nsyn");
	my %matrix = incidence_matrix();
	print_incidence_matrix(\%matrix, "C:/Users/weidewind/workspace/perlCoevolution/TreeUtils/Phylo/MutMap/");
	my $step = 2;
	my $root = $static_tree-> get_root;
	my @array;
	my %hist;
	print "radius,site,node,observed,expected\n";
	my @group;

	foreach my $ind (1..7){
		print "\nindex $ind\n";
		foreach my $node(@{$static_nodes_with_sub{$ind}}){
			
			if(ref($node) eq "REF"){
				$node = ${$node};
			}
			print "\n   node ".$node->get_name()."\n";
			if (!$node->is_terminal){
			my @args = ($ind, $step, $node);
			my_visit_depth_first ($node, \@array,\&update_ring,\&has_no_mutation,\@args,0);
			my_visit_depth_first ($node, \@array,\&max_depth,\&has_no_mutation,\@args,0);
			my $total_muts;
			my $total_length;
	#print "depth ".$static_depth_hash{$ind}{$node->get_name()}."\n";
			if ($static_depth_hash{$ind}{$node->get_name()} > 0){
			foreach my $bin (keys %{$static_ring_hash{$ind}{$node->get_name()}}){
				print "         bin ".$bin*$step."\n";
				print "            adding ".$static_ring_hash{$ind}{$node->get_name()}{$bin}[0]." to muts\n";
				print "            adding ".$static_ring_hash{$ind}{$node->get_name()}{$bin}[1]." to length\n";
				$total_muts += $static_ring_hash{$ind}{$node->get_name()}{$bin}[0];
				$total_length += $static_ring_hash{$ind}{$node->get_name()}{$bin}[1];
			}
			#if ($total_length > 0 && $total_muts/$total_length < 0.005){
			if ($total_length > 0){
			foreach my $bin (sort {$a <=> $b} (keys %{$static_ring_hash{$ind}{$node->get_name()}})){
				
				if ($total_length > 0 && $static_ring_hash{$ind}{$node->get_name()}{$bin}[1] > 0){ #there are some internal nodes with 0-length terminal daughter branches
					print "         bin ".$bin*$step."\n";
					$hist{$bin}[0] += $static_ring_hash{$ind}{$node->get_name()}{$bin}[0]; #observed
					print "            adding ".$static_ring_hash{$ind}{$node->get_name()}{$bin}[0]." to observed\n";
					$hist{$bin}[1] += $total_muts*$static_ring_hash{$ind}{$node->get_name()}{$bin}[1]/$total_length; #expected
					print "            adding ".$total_muts*$static_ring_hash{$ind}{$node->get_name()}{$bin}[1]/$total_length." to expected\n";
				}
				if (!$hist{$bin}[0]){
					$hist{$bin}[0] += 0;
				}
				if (!$hist{$bin}[1]){
					$hist{$bin}[1] += 0;
				}
				#print "$bin,$ind,".$node->get_name().",".$hist{$bin}[0].",".$hist{$bin}[1]."\n";
			}
			}
			}
		}
		}
	}
	
	foreach my $bin (sort {$a <=> $b} keys %hist){
		print " up to ".$bin*$step."\t".$hist{$bin}[0]."\t".$hist{$bin}[1]."\n";
	}
}



sub depth_groups_entrenchment {
	my $step = $_[0];
	my $restriction = $_[1];
	my $root = $static_tree-> get_root;
	my @array;
	my %hist;
	print "radius,site,node,observed,expected\n";
	my @group;
	if ($_[2]){
		@group = @{$_[2]};
	}
	else {
		@group = (1..565);
	}
	foreach my $ind (@group){
		foreach my $node(@{$static_nodes_with_sub{$ind}}){
			if(ref($node) eq "REF"){
				$node = ${$node};
			}
			if (!$node->is_terminal){
			my @args = ($ind, $step, $node);
			my_visit_depth_first ($node, \@array,\&update_ring,\&has_no_mutation,\@args,0);
			my_visit_depth_first ($node, \@array,\&max_depth,\&has_no_mutation,\@args,0);
			my $total_muts;
			my $total_length;
	#print "depth ".$static_depth_hash{$ind}{$node->get_name()}."\n";
			if ($static_depth_hash{$ind}{$node->get_name()} > $restriction){
			foreach my $bin (keys %{$static_ring_hash{$ind}{$node->get_name()}}){
				$total_muts += $static_ring_hash{$ind}{$node->get_name()}{$bin}[0];
				$total_length += $static_ring_hash{$ind}{$node->get_name()}{$bin}[1];
			}
			print $node->get_name()." ".$ind." TOTALS: $total_muts, $total_length, maxdepth ".$static_depth_hash{$ind}{$node->get_name()}."\n";
			#if ($total_length > 0 && $total_muts/$total_length < 0.005){
			if ($total_length > 0){
			foreach my $bin (sort {$a <=> $b} (keys %{$static_ring_hash{$ind}{$node->get_name()}})){
				
				if ($total_length > 0 && $static_ring_hash{$ind}{$node->get_name()}{$bin}[1] > 0){ #there are some internal nodes with 0-length terminal daughter branches
					$hist{$bin}[0] += $static_ring_hash{$ind}{$node->get_name()}{$bin}[0]; #observed
					$hist{$bin}[1] += $total_muts*$static_ring_hash{$ind}{$node->get_name()}{$bin}[1]/$total_length; #expected
					
				}
				if (!$hist{$bin}[0]){
					$hist{$bin}[0] += 0;
				}
				if (!$hist{$bin}[1]){
					$hist{$bin}[1] += 0;
				}
				#print "$bin,$ind,".$node->get_name().",".$hist{$bin}[0].",".$hist{$bin}[1]."\n";
			}
			}
			}
		}
		}
	}
	
	foreach my $bin (sort {$a <=> $b} keys %hist){
		print " up to ".$bin*$step."\t".$hist{$bin}[0]."\t".$hist{$bin}[1]."\n";
	}
	return %hist;
	
	
}

sub no_check{
	return 1;
}


sub enrichment_blue_violet_optimized {
	my $obsfile = $_[0];
	my $bootfile = $_[1];
	my @group = @{$_[2]};
	
	my %grouphash;
	foreach my $ind(@group){
		$grouphash{$ind} = 1;
	}
	
	my %boot_obshash;
	my %boot_exphash;
	my %obshash;
	my %exphash;

	
	open FILE, "<$obsfile" or die "Cannot open $obsfile\n";
	
	my $maxbin = 0;
	while(<FILE>){
		if ($_ =~ /^rad/){
			next;
		}
		my @arr = split(/,/);
		my $ind = $arr[1] ;
		my $bin = $arr[0];
		
		$maxbin = max($bin, $maxbin);
		# bin, ind, node, obs, exp
		$obshash{$ind}{$bin} += $arr[3];
		$exphash{$ind}{$bin} += $arr[4];
	}

	my $norm;
	foreach my $bin(1..$maxbin){
		foreach my $ind (keys %obshash){
			$norm += $obshash{$ind}{$bin};
		}
	}
	
	my %obshash_g;
	my %obshash_c;
	my %exphash_g;
	my %exphash_c;
	
	my $counter_g;
	my $counter_c;
	
	foreach my $bin(1..$maxbin){
			foreach my $ind(keys %obshash){
				if ($grouphash{$ind}){
					$obshash_g{$bin} += $obshash{$ind}{$bin}; # group
					$exphash_g{$bin} += $exphash{$ind}{$bin}; # group
				}
				else {
					$obshash_c{$bin} += $obshash{$ind}{$bin}; # complement
					$exphash_c{$bin} += $exphash{$ind}{$bin}; # complement
				}
			}
	}
		
		
		foreach my $ind(keys %obshash){
				if ($grouphash{$ind}){
					if ($obshash{$ind}){
						$counter_g ++;
					}
				}
				else {
					if ($obshash{$ind}){
						$counter_c ++;
					}
				}
			}
		
		
	print "number of sites: group $counter_g , complement $counter_c\n";
		
	my $obs_median_g = hist_median_for_hash(\%obshash_g);
	my $obs_median_c = hist_median_for_hash(\%obshash_c);
	my $exp_median_g = hist_median_for_hash(\%exphash_g);
	my $exp_median_c = hist_median_for_hash(\%exphash_c);
	print " obs_median_g $obs_median_g  obs_median_c $obs_median_c exp_median_g $exp_median_g  exp_median_c $exp_median_c\n ";	
	
	open FILE, "<$bootfile" or die "Cannot open $bootfile\n";
	my $iteration = 0;
	while(<FILE>){
		if ($_ =~ /^rad/){
			$iteration++;
			next;
		}
		my @arr = split(/,/);
		my $ind = $arr[1] ;
		my $bin = $arr[0];
		# bin, ind, node, obs,exp
		$boot_obshash{$iteration}{$ind}{$bin} += $arr[3];
		$boot_exphash{$iteration}{$ind}{$bin} += $arr[4];
	}
	close FILE;
	print "$iteration iterations\n";
	my $pval_env_enrichment;
	my $pval_env_depletion;
	my $pval_epi_enrichment;
	my $pval_epi_depletion;
	

	
	
	foreach my $i(keys %boot_obshash){
		my $sum;
		my %boot_obshash_g;
		my %boot_obshash_c;
		my %boot_exphash_g;
		my %boot_exphash_c;
		
		my $counter_g;
		my $counter_c;
			foreach my $ind(keys %{ $boot_obshash{$i}}){
				if ($grouphash{$ind}){
					if ($boot_obshash{$i}{$ind}){
						$counter_g ++;
					}
				}
				else {
					if ($boot_obshash{$i}{$ind}){
						$counter_c ++;
					}
				}
			}
		
		
	print "number of sites: group $counter_g , complement $counter_c\n";
		
		
		
		foreach my $bin(1..$maxbin){
			foreach my $ind (1..565){
				$sum += $boot_obshash{$i}{$ind}{$bin};
			}
		}
		foreach my $bin(1..$maxbin){
			foreach my $ind(keys %{$boot_obshash{$i}}){
				if ($grouphash{$ind}){
					$boot_obshash_g{$bin} += $boot_obshash{$i}{$ind}{$bin}*$norm/$sum; # group
					$boot_exphash_g{$bin} += $boot_exphash{$i}{$ind}{$bin}*$norm/$sum; # group
				}
				else {
					$boot_obshash_c{$bin} += $boot_obshash{$i}{$ind}{$bin}*$norm/$sum; # complement
					$boot_exphash_c{$bin} += $boot_exphash{$i}{$ind}{$bin}*$norm/$sum; # complement
				}
			}
		}
		my $boot_obs_median_g = hist_median_for_hash(\%boot_obshash_g);
		my $boot_obs_median_c = hist_median_for_hash(\%boot_obshash_c);
		my $boot_exp_median_g = hist_median_for_hash(\%boot_exphash_g);
		my $boot_exp_median_c = hist_median_for_hash(\%boot_exphash_c);
		print " boot_obs_median_g $boot_obs_median_g  boot_obs_median_c $boot_obs_median_c boot_exp_median_g $boot_exp_median_g  boot_exp_median_c $boot_exp_median_c\n ";	
		
		if (($boot_obs_median_g-$boot_exp_median_g) - ($boot_obs_median_c-$boot_exp_median_c) 
			 >= ($obs_median_g-$exp_median_g) - ($obs_median_c-$exp_median_c) ){
			$pval_env_enrichment += 1;
		}
		if (($boot_obs_median_g-$boot_exp_median_g) - ($boot_obs_median_c-$boot_exp_median_c) 
			 <= ($obs_median_g-$exp_median_g) - ($obs_median_c-$exp_median_c)){
			$pval_env_depletion += 1;
		}
		if (-($boot_obs_median_g-$boot_exp_median_g) + ($boot_obs_median_c-$boot_exp_median_c) 
			 >= -($obs_median_g-$exp_median_g) + ($obs_median_c-$exp_median_c)){
			$pval_epi_enrichment += 1;
		}
		if (-($boot_obs_median_g-$boot_exp_median_g) + ($boot_obs_median_c-$boot_exp_median_c) 
			 <= -($obs_median_g-$exp_median_g) + ($obs_median_c-$exp_median_c)){
			$pval_epi_depletion += 1;
		}
		
	}
	
	print " pval epi enrichment ".$pval_epi_enrichment/$iteration."\n";
	print " pval epi depletion ".$pval_epi_depletion/$iteration."\n";
	print " pval env enrichment ".$pval_env_enrichment/$iteration."\n";
	print " pval env depletion ".$pval_env_depletion/$iteration."\n";
}


sub enrichment_optimized {
	my $obsfile = $_[0];
	my $bootfile = $_[1];
	my @group = @{$_[2]};
	
	my %grouphash;
	foreach my $ind(@group){
		$grouphash{$ind} = 1;
	}
	
	my %hash;
	my %obshash;

	
	open FILE, "<$obsfile" or die "Cannot open $obsfile\n";
	my $maxbin = 0;
	while(<FILE>){
		if ($_ =~ /^rad/){
			next;
		}
		my @arr = split(/,/);
		my $ind = $arr[1] ;
		my $bin = $arr[0];
		$maxbin = max($bin, $maxbin);
		# bin, ind, node, obs
		$obshash{$ind}{$bin} += $arr[3];
	}

	my $norm;
	foreach my $bin(1..$maxbin){
		foreach my $ind (keys %obshash){
			$norm += $obshash{$ind}{$bin};
		}
	}
	
	my %obshash_g;
	my %obshash_c;
	foreach my $bin(1..$maxbin){
			foreach my $ind(keys %obshash){
				if ($grouphash{$ind}){
					$obshash_g{$bin} += $obshash{$ind}{$bin}; # group
				}
				else {
					$obshash_c{$bin} += $obshash{$ind}{$bin}; # complement
				}
			}
	}
		
	my $obs_median_g = hist_median_for_hash(\%obshash_g);
	my $obs_median_c = hist_median_for_hash(\%obshash_c);
	print " obs_median_g $obs_median_g  obs_median_c $obs_median_c\n ";	
	
	open FILE, "<$bootfile" or die "Cannot open $bootfile\n";
	my $iteration = 0;
	while(<FILE>){
		if ($_ =~ /^rad/){
			$iteration++;
			next;
		}
		my @arr = split(/,/);
		my $ind = $arr[1] ;
		my $bin = $arr[0];
		# bin, ind, node, obs
		$hash{$iteration}{$ind}{$bin} += $arr[3];
	}
	close FILE;
	print "$iteration iterations\n";
	my $pval_env_enrichment;
	my $pval_env_depletion;
	my $pval_epi_enrichment;
	my $pval_epi_depletion;
	
	
	foreach my $i(keys %hash){
		my $sum;
		my %hash_g;
		my %hash_c;
		foreach my $bin(1..$maxbin){
			foreach my $ind (1..565){
				$sum += $hash{$i}{$ind}{$bin};
			}
		}
		foreach my $bin(1..$maxbin){
			foreach my $ind(keys %{$hash{$i}}){
				if ($grouphash{$ind}){
					$hash_g{$bin} += $hash{$i}{$ind}{$bin}*$norm/$sum; # group
				}
				else {
					$hash_c{$bin} += $hash{$i}{$ind}{$bin}*$norm/$sum; # complement
				}
			}
		}
		my $boot_median_g = hist_median_for_hash(\%hash_g);
		my $boot_median_c = hist_median_for_hash(\%hash_c);
		print " boot_median_g $boot_median_g  boot_median_c $boot_median_c\n ";	
		
		if ($obs_median_g-$obs_median_c <= $boot_median_g-$boot_median_c){
			$pval_env_enrichment += 1;
		}
		if ($obs_median_g-$obs_median_c >= $boot_median_g-$boot_median_c){
			$pval_env_depletion += 1;
		}
		if (-$obs_median_g+$obs_median_c <= -$boot_median_g+$boot_median_c){
			$pval_epi_enrichment += 1;
		}
		if (-$obs_median_g+$obs_median_c >= -$boot_median_g+$boot_median_c){
			$pval_epi_depletion += 1;
		}
		
	}
	
	print " pval epi enrichment ".$pval_epi_enrichment/$iteration."\n";
	print " pval epi depletion ".$pval_epi_depletion/$iteration."\n";
	print " pval env enrichment ".$pval_env_enrichment/$iteration."\n";
	print " pval env depletion ".$pval_env_depletion/$iteration."\n";
		
}

sub depth_groups_entrenchment_optimized {
	my $step = $_[0];
	my $restriction = $_[1];
	my $root = $static_tree-> get_root;
	my @array;
	my %hist;
	print "radius,site,node,observed,expected\n";
	print 
	my @group;
	if ($_[2]){
		@group = @{$_[2]};
	}
	else {
		@group = (1..565);
	}

	my %closest_ancestors;
	$root->set_generic("-closest_ancestors" => \%closest_ancestors);
	my @args = ( $step, $root);
	my_visit_depth_first ($root, \@array,\&entrenchment_visitor,\&no_check,\@args,0);
	foreach my $ind (@group){
		foreach my $node(@{$static_nodes_with_sub{$ind}}){
			if(ref($node) eq "REF"){
				$node = ${$node};
			}
			my $total_muts;
			my $total_length;
	#print "depth ".$static_depth_hash{$ind}{$node->get_name()}."\n";
			if ($static_subtree_info{$node->get_name()}{$ind}{"maxdepth"} > $restriction){
				foreach my $bin (keys %{$static_subtree_info{$node->get_name()}{$ind}{"hash"}}){
					$total_muts += $static_subtree_info{$node->get_name()}{$ind}{"hash"}{$bin}[0];
					$total_length += $static_subtree_info{$node->get_name()}{$ind}{"hash"}{$bin}[1];
				}	
				#print $node->get_name()." ".$ind." TOTALS: $total_muts, $total_length, maxdepth ".$static_subtree_info{$node->get_name()}{$ind}{"maxdepth"}."\n";
				#if ($total_length > 0 && $total_muts/$total_length < 0.005){
				if ($total_length > 0){
					foreach my $bin (sort {$a <=> $b} (keys %{$static_subtree_info{$node->get_name()}{$ind}{"hash"}})){
						if ($total_length > 0 && $static_subtree_info{$node->get_name()}{$ind}{"hash"}{$bin}[1] > 0){ #there are some internal nodes with 0-length terminal daughter branches
							$hist{$bin}[0] += $static_subtree_info{$node->get_name()}{$ind}{"hash"}{$bin}[0]; #observed
							$hist{$bin}[1] += $total_muts*$static_subtree_info{$node->get_name()}{$ind}{"hash"}{$bin}[1]/$total_length; #expected	
						}
						if (!$hist{$bin}[0]){
							$hist{$bin}[0] += 0;
						}
						if (!$hist{$bin}[1]){
							$hist{$bin}[1] += 0;
						}
	 print "$bin,$ind,".$node->get_name().",".$hist{$bin}[0].",".$hist{$bin}[1]."\n";
				}
				}
			}
			
		}
	}	
	
	#foreach my $bin (sort {$a <=> $b} keys %hist){
	#	print " up to ".$bin*$step."\t".$hist{$bin}[0]."\t".$hist{$bin}[1]."\n";
	#}
	return %hist;	
}





sub check_depth_groups_entrenchment_optimized {
	my $step = $_[0];
	my $subtract_tallest = $_[1];
	my $root = $static_tree-> get_root;
	my @array;
	my %hist;
	#print "mutnum,maxdepth,radius,site,node,observed,expected\n";
	print 
	my @group;
	if ($_[2]){
		@group = @{$_[2]};
	}
	else {
		@group = (1..565);
	}

	my %closest_ancestors;
	$root->set_generic("-closest_ancestors" => \%closest_ancestors);
	my @args = ( $step, $root, $subtract_tallest);
	my_visit_depth_first ($root, \@array,\&entrenchment_visitor,\&no_check,\@args,0);
	foreach my $ind (@group){
		foreach my $node(@{$static_nodes_with_sub{$ind}}){
			if(ref($node) eq "REF"){
				$node = ${$node};
			}
			my $total_muts;
			my $total_length;

				#my $maxdepth = int($static_subtree_info{$node->get_name()}{$ind}{"maxdepth"}/50);
				if ($static_subtree_info{$node->get_name()}{$ind}{"maxdepth"} > 0){
				my %subtract_hash;
				
				if ($subtract_tallest){
					my $tallest_tip = ${$static_hash_of_nodes{$static_subtree_info{$node->get_name()}{$ind}{"maxdepth_node"}}};	
					my @path_to_tallest_tip = @{$tallest_tip->get_ancestors()};	
					my $index_of_ancestor_node;
					my $path_length = $tallest_tip->get_branch_length; # todo: check if ancestors contain the node itself
				
					for(my $n = 0; $n < scalar @path_to_tallest_tip; $n++){		
						if ($path_to_tallest_tip[$n] eq $node){
							$index_of_ancestor_node = $n;
							last;
						}
						my $depth = $static_distance_hash{$node->get_name()}{$path_to_tallest_tip[$n]->get_name()};
						$subtract_hash{bin($depth,$step)} += $path_to_tallest_tip[$n]->get_branch_length;
						$path_length += $path_to_tallest_tip[$n]->get_branch_length;
					}
				}
				
					
					
				foreach my $bin (keys %{$static_subtree_info{$node->get_name()}{$ind}{"hash"}}){
					$total_muts += $static_subtree_info{$node->get_name()}{$ind}{"hash"}{$bin}[0];
					$total_length += $static_subtree_info{$node->get_name()}{$ind}{"hash"}{$bin}[1];
					if ($subtract_tallest && $subtract_hash{$bin}){
						$total_length -= $subtract_hash{$bin};
					}
				}	

				if ($total_length > 0){
					$hist{$total_muts}{"counter"} += 1;
					foreach my $bin (sort {$a <=> $b} (keys %{$static_subtree_info{$node->get_name()}{$ind}{"hash"}})){
						if ($total_length > 0 && $static_subtree_info{$node->get_name()}{$ind}{"hash"}{$bin}[1] > 0){ #there are some internal nodes with 0-length terminal daughter branches
							my $local_length = $static_subtree_info{$node->get_name()}{$ind}{"hash"}{$bin}[1];
							if ($subtract_tallest && $subtract_hash{$bin}){
								$local_length -= $subtract_hash{$bin};
							}
							$hist{$total_muts}{$bin}[0] += $static_subtree_info{$node->get_name()}{$ind}{"hash"}{$bin}[0]; #observed
							$hist{$total_muts}{$bin}[1] += $total_muts*$local_length/$total_length; #expected	
							
						}
						if (!$hist{$total_muts}{$bin}[0]){
							$hist{$total_muts}{$bin}[0] = 0;
						}
						if (!$hist{$total_muts}{$bin}[1]){
							$hist{$total_muts}{$bin}[1] = 0;
						}
	# print "$total_muts,".$static_subtree_info{$node->get_name()}{$ind}{"maxdepth"}."$bin,$ind,".$node->get_name().",".$hist{$bin}[0].",".$hist{$bin}[1]."\n";
				}
				}
				
		}
		}
	}	
	
	#foreach my $bin (sort {$a <=> $b} keys %hist){
	#	print " up to ".$bin*$step."\t".$hist{$bin}[0]."\t".$hist{$bin}[1]."\n";
	#}
	return %hist;	
}



sub depth_groups_entrenchment_hash {
	my $step = $_[0];
	my $root = $static_tree-> get_root;
	my @array;
	my %hist;
	print "radius,site,node,observed,expected\n";
	my @group;
	if ($_[1]){
		@group = @{$_[1]};
	}
	else {
		@group = (1..565);
	}
	foreach my $ind (@group){
		foreach my $node(@{$static_nodes_with_sub{$ind}}){
			if(ref($node) eq "REF"){
				$node = ${$node};
			}
			if (!$node->is_terminal){
			my @args = ($ind, $step, $node);
			my_visit_depth_first ($node, \@array,\&update_ring,\&has_no_mutation,\@args,0);
			my_visit_depth_first ($node, \@array,\&max_depth,\&has_no_mutation,\@args,0);
			my $total_muts;
			my $total_length;
	#print "depth ".$static_depth_hash{$ind}{$node->get_name()}."\n";
			if ($static_depth_hash{$ind}{$node->get_name()} > 50){
			foreach my $bin (keys %{$static_ring_hash{$ind}{$node->get_name()}}){
				$total_muts += $static_ring_hash{$ind}{$node->get_name()}{$bin}[0];
				$total_length += $static_ring_hash{$ind}{$node->get_name()}{$bin}[1];
			}
			#if ($total_length > 0 && $total_muts/$total_length < 0.005){
			if ($total_length > 0){
			foreach my $bin (sort {$a <=> $b} (keys %{$static_ring_hash{$ind}{$node->get_name()}})){
				
				if ($total_length > 0 && $static_ring_hash{$ind}{$node->get_name()}{$bin}[1] > 0){ #there are some internal nodes with 0-length terminal daughter branches
					$hist{$bin}{$ind}[0] += $static_ring_hash{$ind}{$node->get_name()}{$bin}[0]; #observed
					$hist{$bin}{$ind}[1] += $total_muts*$static_ring_hash{$ind}{$node->get_name()}{$bin}[1]/$total_length; #expected
					
				}
				if (!$hist{$bin}{$ind}[0]){
					$hist{$bin}{$ind}[0] += 0;
				}
				if (!$hist{$bin}{$ind}[1]){
					$hist{$bin}{$ind}[1] += 0;
				}
				
			}
			}
			}
		}
		}

	}
	
	foreach my $bin (sort {$a <=> $b} keys %hist){
		my $obs_sum;
		my $exp_sum;
		foreach my $s (@group){
			$obs_sum += $hist{$bin}{$s}[0];
			$exp_sum += $hist{$bin}{$s}[1];
		}
		print " up to ".$bin*$step."\t".$obs_sum."\t".$exp_sum."\n";
	}
	return %hist;
	
}


sub depth_groups_entrenchment_bootstrap {
	my $step = $_[0];
	my $root = $static_tree-> get_root;
	my @pregroup = @{$_[1]};
	my $iterations = $_[2];
	
	# both group and complement must contain only variable sites:
	my %variable_sites;
	foreach my $s(1..565){
		if ($static_nodes_with_sub{$s}){
			$variable_sites{$s} = 1;
		}
	}
	
	my %ghash;
	my @group;
	foreach my $s(@pregroup){
		if ($variable_sites{$s}){
			$ghash{$s} = 1;
			push @group, $s;
		}
	}
	my @complement;
	foreach my $s(keys %variable_sites){
		if (!$ghash{$s}){
			push @complement, $s;
		}
	}
	my @variable_sites_array = keys %variable_sites;
	my %hist = depth_groups_entrenchment_hash($step, \@variable_sites_array);
	my @obs_values = diffdiff(\%hist, \@group, \@complement);
	 
	my $counter1;
	my $counter2;

	
	for (my $i = 0; $i < $iterations; $i++){
		
		my @boot_group = shuffle @variable_sites_array;
		my @boot_complement = splice (@boot_group, scalar @group, scalar @variable_sites_array - scalar @group);
		my @boot_values = diffdiff(\%hist, \@boot_group, \@boot_complement);
		
		if ($boot_values[4] >= $obs_values[4]){
			$counter1++;
		}
		if ($boot_values[4] <= $obs_values[4]){
			$counter2++;
		}


	}
	print "epistasis enrichment ".$counter1/$iterations."\n";
	print "environment enrichment ".$counter2/$iterations."\n";

}


sub squashed_hist {
	my %hist = %{$_[0]};
	my %squashed_hist;
	foreach my $bin (keys %hist){
		foreach my $ind (keys %{$hist{$bin}}){
			$squashed_hist{$bin}[0] += $hist{$bin}{$ind}[0];
			$squashed_hist{$bin}[1] += $hist{$bin}{$ind}[1];
		}
	}
	return %squashed_hist;
}

sub squashed_hist_subset {
	my %hist = %{$_[0]};
	my @group = @{$_[1]};
	my %hist_subset;
	foreach my $bin (keys %hist){
		foreach my $ind (@group){
			$hist_subset{$bin}[0] += $hist{$bin}{$ind}[0];
			$hist_subset{$bin}[1] += $hist{$bin}{$ind}[1];
		}
	}
	return %hist_subset;
}

sub diffdiff {
	my %hist = %{$_[0]};
	my @group = @{$_[1]};
	my @complement = @{$_[2]};
	
	my %group_hist = squashed_hist_subset(\%hist, \@group);
	my %group_obs_hist;
	my %group_exp_hist;
	foreach my $k(keys %group_hist){
		$group_obs_hist{$k} = $group_hist{$k}->[0];
		$group_exp_hist{$k} = $group_hist{$k}->[1];
	}
	my %complement_hist = squashed_hist_subset(\%hist, \@complement);
	my %complement_obs_hist;
	my %complement_exp_hist;
	foreach my $k(keys %complement_hist){
		$complement_obs_hist{$k} = $complement_hist{$k}->[0];
		$complement_exp_hist{$k} = $complement_hist{$k}->[1];
	}
	my $group_obs = hist_median_for_hash(\%group_obs_hist);
	my $group_exp = hist_median_for_hash(\%group_exp_hist);
	my $compl_obs = hist_median_for_hash(\%complement_obs_hist);
	my $compl_exp = hist_median_for_hash(\%complement_exp_hist);
	
	print "\n Medians:\t".$group_obs."\t";
	print $group_exp."\t";
	print $compl_obs."\t";
	print $compl_exp."\t";

	my $diffdiff =  $group_exp - $group_obs - $compl_exp + $compl_obs; # > 0 if epist
	print $diffdiff."\n";
	
	return ($group_obs, $group_exp, $compl_obs, $compl_exp, $diffdiff);
}

# incorrect! updates ring (adds to it) with every iteration
sub _delete_depth_groups_entrenchment_bootstrap {
	
	my $step = $_[0];
	my $root = $static_tree-> get_root;
	my @pregroup = @{$_[1]};
	my $iterations = $_[2];
	
	# both group and complement must contain only variable sites:
	my %variable_sites;
	foreach my $s(1..565){
		if ($static_nodes_with_sub{$s}){
			$variable_sites{$s} = 1;
		}
	}
	
	my %ghash;
	my @group;
	foreach my $s(@pregroup){
		if ($variable_sites{$s}){
			$ghash{$s} = 1;
			push @group, $s;
		}
	}
	my @complement;
	foreach my $s(keys %variable_sites){
		if (!$ghash{$s}){
			push @complement, $s;
		}
	}
	
	my %group_hist = depth_groups_entrenchment($step, \@group);
	my %group_obs_hist;
	my %group_exp_hist;
	foreach my $k(keys %group_hist){
		$group_obs_hist{$k} = $group_hist{$k}->[0];
		$group_exp_hist{$k} = $group_hist{$k}->[1];
	}
	my %complement_hist = depth_groups_entrenchment($step, \@complement);
	my %complement_obs_hist;
	my %complement_exp_hist;
	foreach my $k(keys %complement_hist){
		$complement_obs_hist{$k} = $complement_hist{$k}->[0];
		$complement_exp_hist{$k} = $complement_hist{$k}->[1];
	}
	my $group_obs = hist_median_for_hash(\%group_obs_hist);
	my $group_exp = hist_median_for_hash(\%group_exp_hist);
	my $compl_obs = hist_median_for_hash(\%complement_obs_hist);
	my $compl_exp = hist_median_for_hash(\%complement_exp_hist);
	
	
	print "\n Medians:\t".$group_obs."\t";
	print $group_exp."\t";
	print $compl_obs."\t";
	print $compl_exp."\t";

	my $diffdiff =  $group_exp - $group_obs - $compl_exp + $compl_obs; # > 0 if epist
	 
	my @variable_sites_array = keys %variable_sites;
	my $counter1;
	my $counter2;
	my $counter3;
	my $counter4;
	
	for (my $i = 0; $i < $iterations; $i++){
		
		my @boot_group = shuffle @variable_sites_array;
		my @boot_complement = splice (@boot_group, scalar @group, scalar @variable_sites_array - scalar @group);
		foreach my $g(@boot_group){
			print $g."\n";
		}
		
		
		my %boot_group_hist = depth_groups_entrenchment($step, \@boot_group);
		my %boot_group_obs_hist;
		my %boot_group_exp_hist;
		foreach my $k(keys %boot_group_hist){
			$boot_group_obs_hist{$k} = $boot_group_hist{$k}->[0];
			$boot_group_exp_hist{$k} = $boot_group_hist{$k}->[1];
		}
		my %boot_complement_hist = depth_groups_entrenchment($step, \@boot_complement);
		my %boot_complement_obs_hist;
		my %boot_complement_exp_hist;
		foreach my $k(keys %boot_complement_hist){
			$boot_complement_obs_hist{$k} = $boot_complement_hist{$k}->[0];
			$boot_complement_exp_hist{$k} = $boot_complement_hist{$k}->[1];
		}
		my $boot_group_obs = hist_median_for_hash(\%boot_group_obs_hist);
		my $boot_group_exp = hist_median_for_hash(\%boot_group_exp_hist);
		my $boot_compl_obs = hist_median_for_hash(\%boot_complement_obs_hist);
		my $boot_compl_exp = hist_median_for_hash(\%boot_complement_exp_hist);
		
		if ($boot_group_exp-$boot_group_obs - $boot_compl_exp+$boot_compl_obs >= $diffdiff){
			$counter1++;
		}
		if ($boot_group_exp-$boot_group_obs - $boot_compl_exp+$boot_compl_obs <= $diffdiff){
			$counter2++;
		}
		if (-$boot_group_exp+$boot_group_obs + $boot_compl_exp-$boot_compl_obs >= $diffdiff){
			$counter3++;
		}
		if (-$boot_group_exp+$boot_group_obs + $boot_compl_exp-$boot_compl_obs <= $diffdiff){
			$counter4++;
		}

	}
	print "epistasis enrichment ".$counter1/$iterations."\n";
	print "epistasis depletion ".$counter2/$iterations."\n";
	print "environment enrichment ".$counter3/$iterations."\n";
	print "environment depletion ".$counter4/$iterations."\n";
}

sub sites_num {
	my $counter;
	foreach my $ind (1..566){
		if (ref($static_nodes_with_sub{$ind}) eq "ARRAY") {
		my $c = scalar @{$static_nodes_with_sub{$ind}};
		if ($c > 2){
			$counter++;
		}
		}
	}
	print $counter;
}


sub depth_groups_entrenchment_heaps {
	my $step = $_[0];
	my $root = $static_tree-> get_root;
	my @array;
	my %hist;
	print "radius,site,node,observed,expected\n";
	my @group;
	if ($_[1]){
		@group = @{$_[1]};
	}
	else {
		@group = (1..565);
	}
	foreach my $ind (@group){
		foreach my $node(@{$static_nodes_with_sub{$ind}}){
			$node = ${$node};
			if (!$node->is_terminal){
			my @args = ($ind, $step, $node);
			my_visit_depth_first ($node, \@array,\&update_ring,\&has_no_mutation,\@args,0);
			my_visit_depth_first ($node, \@array,\&max_depth,\&has_no_mutation,\@args,0);
			my $total_muts;
			my $total_length;
	#print "depth ".$static_depth_hash{$ind}{$node->get_name()}."\n";
			if ($static_depth_hash{$ind}{$node->get_name()} > 50){
			foreach my $bin (keys %{$static_ring_hash{$ind}{$node->get_name()}}){
				$total_muts += $static_ring_hash{$ind}{$node->get_name()}{$bin}[0];
				$total_length += $static_ring_hash{$ind}{$node->get_name()}{$bin}[1];
			}
			if ($total_length > 0 && $total_muts/$total_length < 0.005){
			#if ($total_length > 0){
			foreach my $bin (sort {$a <=> $b} (keys %{$static_ring_hash{$ind}{$node->get_name()}})){
				my $heap0 = 0;
				my $heap1 = 1;
				if ($static_ring_hash{$ind}{$node->get_name()}{1}[0] > $total_muts*$static_ring_hash{$ind}{$node->get_name()}{1}[1]/$total_length) {
					$heap0 = 2;
					$heap1 = 3;
				}
				if ($total_length > 0 && $static_ring_hash{$ind}{$node->get_name()}{$bin}[1] > 0){ #there are some internal nodes with 0-length terminal daughter branches
					$hist{$bin}[$heap0] += $static_ring_hash{$ind}{$node->get_name()}{$bin}[0]; #observed
					$hist{$bin}[$heap1] += $total_muts*$static_ring_hash{$ind}{$node->get_name()}{$bin}[1]/$total_length; #expected
					
				}
				if (!$hist{$bin}[$heap0]){
					$hist{$bin}[$heap0] += 0;
				}
				if (!$hist{$bin}[$heap1]){
					$hist{$bin}[$heap1] += 0;
				}
				print "$bin,$ind,".$node->get_name().",".$hist{$bin}[0].",".$hist{$bin}[1]."\n";
			}
			}
			}
		}
		}
	}
	
	foreach my $bin (sort {$a <=> $b} keys %hist){
		print " up to ".$bin*$step."\t".$hist{$bin}[0]."\t".$hist{$bin}[1]."\t".$hist{$bin}[2]."\t".$hist{$bin}[3]."\n";
	}
	
	
}


# one hist for one ancestor aa in one site
sub egor_site_entrenchment {
	my $step = $_[0];
	my $root = $static_tree-> get_root;
	my @array;
	print "radius,site,node,density\n";
	for (my $ind = 1; $ind < 566; $ind++){
		foreach my $node(@{$static_nodes_with_sub{$ind}}){
			$node = ${$node};
			if (!$node->is_terminal){
			my %hist;
			#my $ancestor = ${$static_subs_on_node{$node->get_name()}}{$ind}->{"Substitution::ancestral_allele"};
			my @args = ($ind, $step, $node);
			my_visit_depth_first ($node, \@array,\&update_ring,\&has_no_mutation,\@args,0);
			
			my $cumulative_muts;
			my $cumulative_length;
			foreach my $bin (sort {$a <=> $b} (keys %{$static_ring_hash{$ind}{$node->get_name()}})){
				#print "bin $bin observed ".$static_ring_hash{$ind}{$node->get_name()}{$bin}[0]." totmut $total_muts totlen $total_length\n";
				my $muts_in_bin = $static_ring_hash{$ind}{$node->get_name()}{$bin}[0];
				my $length_of_bin = $static_ring_hash{$ind}{$node->get_name()}{$bin}[1];
				if ($length_of_bin > 0){ #there are some internal nodes with 0-length terminal daughter branches
					#	$hist{$bin}{$ancestor} += $static_ring_hash{$ind}{$node->get_name()}{$bin}[0]/$static_ring_hash{$ind}{$node->get_name()}{$bin}[1]; #density
					$cumulative_muts += $muts_in_bin;
					$cumulative_length += $length_of_bin;
					$hist{$bin} = $cumulative_muts/$cumulative_length; #density
				}
				else {
					$hist{$bin} = 0;
				}
				print "$bin,$ind,".$node->get_name().",".$hist{$bin}."\n";
			}
			

		}
		}
	}
	
	
	
	
}

# one hist for one ancestor aa in one site
sub egor_smart_site_entrenchment {
	my $step = 1;
	my $root = $static_tree-> get_root;
	my @array;
	print "radius,site,node,density\n";
	for (my $ind = 1; $ind < 566; $ind++){
		foreach my $node(@{$static_nodes_with_sub{$ind}}){
			$node = ${$node};
			if (!$node->is_terminal){
			my %hist;
			#my $ancestor = ${$static_subs_on_node{$node->get_name()}}{$ind}->{"Substitution::ancestral_allele"};
			my @args = ($ind, $step, $node);
			my_visit_depth_first ($node, \@array,\&update_ring,\&has_no_mutation,\@args,0);
			
			my $cumulative_muts;
			my $cumulative_length;
			my @sorted_keys = sort {$a <=> $b} (keys %{$static_ring_hash{$ind}{$node->get_name()}});
			foreach my $bin (@sorted_keys){
				#print "bin $bin observed ".$static_ring_hash{$ind}{$node->get_name()}{$bin}[0]." totmut $total_muts totlen $total_length\n";
				my $muts_in_bin = $static_ring_hash{$ind}{$node->get_name()}{$bin}[0];
				my $length_of_bin = $static_ring_hash{$ind}{$node->get_name()}{$bin}[1];
				$cumulative_muts += $muts_in_bin;
				$cumulative_length += $length_of_bin;
				if ($cumulative_length > 0){ #there are some internal nodes with 0-length terminal daughter branches
					#	$hist{$bin}{$ancestor} += $static_ring_hash{$ind}{$node->get_name()}{$bin}[0]/$static_ring_hash{$ind}{$node->get_name()}{$bin}[1]; #density

					$hist{$bin} = $cumulative_muts/$cumulative_length; #density
				}
				else {
					$hist{$bin} = 0;
				}
				#if ($muts_in_bin > 0 || $bin == $sorted_keys[0] || $bin == $sorted_keys[-1]){	
					if ($muts_in_bin > 0){	
					print "$bin,$ind,".$node->get_name().",".$hist{$bin}."\n";
				}
				
			}
			

		}
		}
	}
	
	
	
	
}

# last egor plots sub
sub egor_diff_rings_site_entrenchment {
	my $step = 1;
	my $root = $static_tree-> get_root;
	my @array;
	print "radius,site,node,density,cum_muts,cum_length\n";
	for (my $ind = 1; $ind < 566; $ind++){
		foreach my $node(@{$static_nodes_with_sub{$ind}}){
			$node = ${$node};
			if (!$node->is_terminal){
			my %hist;
			#my $ancestor = ${$static_subs_on_node{$node->get_name()}}{$ind}->{"Substitution::ancestral_allele"};
			my @args = ($ind, $step, $node);
			my_visit_depth_first ($node, \@array,\&update_ring,\&has_no_mutation,\@args,0);
			my $cumulative_muts;
			my $cumulative_length;
			my @sorted_keys = sort {$a <=> $b} (keys %{$static_ring_hash{$ind}{$node->get_name()}});
			foreach my $bin (@sorted_keys){
				#print "bin $bin observed ".$static_ring_hash{$ind}{$node->get_name()}{$bin}[0]." totmut $total_muts totlen $total_length\n";
				my $muts_in_bin = $static_ring_hash{$ind}{$node->get_name()}{$bin}[0];
				my $length_of_bin = $static_ring_hash{$ind}{$node->get_name()}{$bin}[1];
				$cumulative_muts += $muts_in_bin;
				$cumulative_length += $length_of_bin;
				if ($cumulative_length > 0){ #there are some internal nodes with 0-length terminal daughter branches
					#	$hist{$bin}{$ancestor} += $static_ring_hash{$ind}{$node->get_name()}{$bin}[0]/$static_ring_hash{$ind}{$node->get_name()}{$bin}[1]; #density
					$hist{$bin} = $cumulative_muts/$cumulative_length; #density
				}
				else {
					$hist{$bin} = 0;
				}
				#if ($muts_in_bin > 0 || $bin == $sorted_keys[0] || $bin == $sorted_keys[-1]){	
					if ($muts_in_bin > 0){	
					print "$bin,$ind,".$node->get_name().",".$hist{$bin}.",".$cumulative_muts.",".$cumulative_length."\n";
					$cumulative_muts = 0;
					$cumulative_length = 0;
				}
				
			}
			

		}
		}
	}
	
	
	
	
}

sub egor_rings_site_entrenchment {
		my $step = 10;
	my $root = $static_tree-> get_root;
	my @array;
	print "radius,site,node,density\n";
	for (my $ind = 1; $ind < 566; $ind++){
		foreach my $node(@{$static_nodes_with_sub{$ind}}){
			$node = ${$node};
			if (!$node->is_terminal){
			my %hist;
			#my $ancestor = ${$static_subs_on_node{$node->get_name()}}{$ind}->{"Substitution::ancestral_allele"};
			my @args = ($ind, $step, $node);
			my_visit_depth_first ($node, \@array,\&update_ring,\&has_no_mutation,\@args,0);

			my @sorted_keys = sort {$a <=> $b} (keys %{$static_ring_hash{$ind}{$node->get_name()}});
			foreach my $bin (@sorted_keys){
				#print "bin $bin observed ".$static_ring_hash{$ind}{$node->get_name()}{$bin}[0]." totmut $total_muts totlen $total_length\n";
				my $muts_in_bin = $static_ring_hash{$ind}{$node->get_name()}{$bin}[0];
				my $length_of_bin = $static_ring_hash{$ind}{$node->get_name()}{$bin}[1];
				if ($length_of_bin > 0){ #there are some internal nodes with 0-length terminal daughter branches
					#	$hist{$bin}{$ancestor} += $static_ring_hash{$ind}{$node->get_name()}{$bin}[0]/$static_ring_hash{$ind}{$node->get_name()}{$bin}[1]; #density
					$hist{$bin} = $muts_in_bin/$length_of_bin; #density
				}
				else {
					$hist{$bin} = 0;
				}
				#if ($muts_in_bin > 0 || $bin == $sorted_keys[0] || $bin == $sorted_keys[-1]){	
				#	if ($muts_in_bin > 0){	
					print "$bin,$ind,".$node->get_name().",".$hist{$bin}."\n";
				#}
				
			}
			

		}
		}
	}
	
	
	
	
}

# input: name of the protein (for filepath), $tree.
# returns a hash: key - node of the tree, value - total length of all branches of its subtree
sub subtree_lengths {
	my $protein_name = $_[0];
	my $tree = $_[1];
	my $file = "C:/Users/weidewind/Documents/CMD/Coevolution/Influenza/".$_[0]."_subtree_hash";
	my $hash;
	if (-e $file){
		$hash = retrieve($file);
	}
	else {
		$tree->visit_depth_first(
			-post_daughter  => sub {	my $node=shift; 
										my $subtree_length = $node->get_first_daughter()->get_generic("-subtree_length") + $node->get_first_daughter()->get_branch_length;
										$subtree_length += $node->get_last_daughter()->get_generic("-subtree_length") + $node->get_last_daughter()->get_branch_length;
										$node->set_generic("-subtree_length" => $subtree_length);
										my $name = $node->get_name;
										$hash->{$name} = $subtree_length;
										print "post_daughter: $name, $subtree_length\n" }
		);
		store (\%{$hash}, $file);
	}
	
	return $hash;
}


sub distances_to_subtree_mutations {
	my $tree = $_[0];
	my $node = $_[1]; # a node, not a nodename
	my %subs_on_node = %{$_[2]};
	my $site_index = $_[3];
	my $subtree = ${$node}->get_subtree;
	$subtree->visit_depth_first(
			-pre_daughter => sub {	my $node=shift; 
									if (${$subs_on_node{${$node}->get_name()}}{$site_index}){
										my $node = shift;
									}
		
			}
			
		);
}

 sub my_visit_depth_first {
		my $node = $_[0];
		my @array = @{$_[1]};
		my $action_callback = $_[2];
		my $check_callback = $_[3];
		my $callback_args = $_[4];
		my $depth = $_[5];
		push @array, $node;
		my $len = $node -> get_branch_length;
		$depth += $len;
		&$action_callback($node, $callback_args, $depth);
		if (! $node->is_terminal && &$check_callback($node, $callback_args)){
			my $i = 0;
			while($node->get_child($i)){
				@array = my_visit_depth_first($node->get_child($i), \@array, \&$action_callback, \&$check_callback, $callback_args, $depth);
				$i++;
				#@array = my_visit_depth_first($node->get_first_daughter, \@array, \&$action_callback, \&$check_callback, $callback_args, $depth);
				#@array = my_visit_depth_first($node->get_last_daughter, \@array, \&$action_callback, \&$check_callback, $callback_args, $depth);
			}
		}
		
		$node = pop @array;
#print $node->get_name()."\t".$depth."\n";
		$depth -= $len;
		return @array;
    }
    
     sub collect_ring_squares {
		my $node = $_[0];
		my $site_index = $_[1]->[0];
		my $step = $_[1]->[1];

		my $root = $static_tree -> get_root;
		my_visit_depth_first($root, (), \&has_no_mutation, \&update_ring, \($site_index, $step));
    }
    
    # old version, does not account for  mutations of the other type
 	sub has_no_same_type_mutation { #has_no_same_type_mutation
 		my $node = $_[0];
 		my $site_index = $_[1]->[0];
 		my $starting_node = $_[1]->[2];
 		
 		if ($node eq $starting_node){
 			return 1;
 		}
 		if (${$static_subs_on_node{$node->get_name()}}{$site_index}){
 			return 0;
 		}
 		else {
 			return 1;
 		}
 	}  
 	
 	
 	# also accounts for mutations of the other type (synonimous for nsyn and non-synonimous for syn)
 	sub has_no_mutation{
 		my $node = $_[0];
 		my $site_index = $_[1]->[0];
 		my $starting_node = $_[1]->[2];
 		
 		if ($node eq $starting_node){
 			return 1;
 		}
 		if (${$static_subs_on_node{$node->get_name()}}{$site_index}){
 			return 0;
 		}
 		
 		if ($static_state eq "n"){
 			if (is_neighbour_changing(${$static_background_subs_on_node{$node->get_name()}}{$site_index}, 1) == 1){
 				return 0;
 			}
 			else {
 				return 1;
 			}
 		}
 		else {
 			if (${$static_background_subs_on_node{$node->get_name()}}{$site_index}){
 				return 0;
 			}
 			else {
 				return 1;
 			}
 		}
 	}  
 	
 	
 	sub has_no_background_mutation {
 		my $node = $_[0];
 		my $site_index = $_[1];
 		if ($static_state eq "n"){
 			if (is_neighbour_changing(${$static_background_subs_on_node{$node->get_name()}}{$site_index}, 1) == 1){
 				return 0;
 			}
 			else {
 				return 1;
 			}
 		}
 		else {
 			if (${$static_background_subs_on_node{$node->get_name()}}{$site_index}){
 				return 0;
 			}
 			else {
 				return 1;
 			}
 		}
 	}
 	
# 	test_max_depth();
 	sub test_max_depth{
 		foreach my $ind (3..3){
 		foreach my $node(@{$static_nodes_with_sub{$ind}}){
			$node = ${$node};
			if (!$node->is_terminal){
				my @array;
				my @args = ($ind, 1, $node);
				print "Now measuring max depth for ".$node->get_name()."\n";
				my_visit_depth_first ($node, \@array,\&max_depth,\&has_no_mutation,\@args,0);
				print $node->get_name()."\t".$static_depth_hash{$ind}{$node->get_name()}."\n";
			}
 		}
 		}
 	}
 	
 	sub max_depth {
 		my $node = $_[0];
		my $site_index = $_[1]->[0];
		my $starting_node = $_[1]->[2];
		my $depth = $_[2] - $starting_node->get_branch_length;
	#	print " current node is ".$node-> get_name()."\n";
	#	print " current depth is $depth\n";
		if ($node eq $starting_node){
	#		print $node-> get_name()." equal ";
			return;
		}
		if ($starting_node -> is_terminal){
	#		print $node-> get_name()." terminal ";
			return;
		}
 		
 			if ($static_depth_hash{$site_index}{$starting_node->get_name()}){
 	#			print "\n".$static_depth_hash{$site_index}{$starting_node->get_name()}.", ".$depth."; max is ";
 				$static_depth_hash{$site_index}{$starting_node->get_name()} = max($static_depth_hash{$site_index}{$starting_node->get_name()}, $depth);
 	#			print $static_depth_hash{$site_index}{$starting_node->get_name()}."\n";
 			}
 			else {
 	#			print "\n init: $depth\n";
 				$static_depth_hash{$site_index}{$starting_node->get_name()} = $depth;
 			}
 		
 	}
 	
 	sub max ($$) { $_[$_[0] < $_[1]] }
 	
 	sub update_ring {
 		my $node = $_[0];
		my $site_index = $_[1]->[0];
		my $step = $_[1]->[1];
		my $starting_node = $_[1]->[2];
		my $depth = $_[2] - $starting_node->get_branch_length ;
		if ($node eq $starting_node){
			#print " \n equality: ".$starting_node ->get_name."\t".$node ->get_name."\n";
			return;
		}
		if ($starting_node -> is_terminal){
			#print " \n terminal: ".$starting_node ->get_name."\n";
			return;
		}
 		#print "depth $depth step $step bin ".(bin($depth,$step))."\n";
 		if (!has_no_same_type_mutation($_[0], \@{$_[1]})){ #has_no_same_type_mutation
 			$static_ring_hash{$site_index}{$starting_node->get_name()}{bin($depth,$step)}[0] += 1;
 		}
 		$static_ring_hash{$site_index}{$starting_node->get_name()}{bin($depth,$step)}[1] += $node->get_branch_length;
 		
 	}
 	
 	#sub set_distance_matrix {
 	#	my $tree = $static_tree;
 	#	my @nodes = @{$tree->get_internals};
 	#	push @nodes, @{$tree->get_terminals};
 	#	foreach my $node(@nodes){
 	#		foreach my $other_node(@nodes){
 	#			$static_distance_hash{$node->get_name()}{$other_node->get_name()} = calc_true_patristic_distance(\${$node}, \${$other_node});
 	#		}
 	#	}
 			
 	#}
 	
 	sub set_dr_distance_matrix {
 		my $prot = $_[0];
 		my $file = "C:/Users/weidewind/Documents/CMD/Coevolution/Influenza/Kryazhimsky11/Mock/".$prot."_distance_matrix.csv";
 		
 		open CSV, "<$file";
 		my $header = <CSV>;
 		$header =~ s/[\s\n\r\t]+$//s;
 		my @nodelables = split(',', $header);
 		#shift @nodelables;
 		while(<CSV>){
 			$_ =~ s/[\s\n\r\t]+$//s;
 			print $_;
 			my @dists = split(',', $_);
 			my $node = $dists[0];
 			for (my $i = 1; $i < scalar @dists; $i++){
 				$static_distance_hash{$node}{$nodelables[$i]} = $dists[$i];
 			#	print " $node ^".$nodelables[$i]."^".$dists[$i]."^";
 			}
 		}
 		close CSV;

 	}
 	
 	
sub set_distance_matrix {
 		my $prot = $_[0];
 		my $file = "C:/Users/weidewind/Documents/CMD/Coevolution/Influenza/Kryazhimsky11/".$prot."_distance_matrix.csv";
 		
 		open CSV, "<$file";
 		my $header = <CSV>;
 		$header =~ s/[\s\n\r\t]+$//s;
 		my @nodelables = split(',', $header);
 		while(<CSV>){
			$_ =~ s/[\s\n\r\t]+$//s;
 			my @dists = split(',', $_);
 			my $node = $dists[0];
 			for (my $i = 1; $i < scalar @dists; $i++){
 				$static_distance_hash{$node}{$nodelables[$i]} = $dists[$i];
 			}
 		}
 		close CSV;

 	}
 	
 	# track_tallest is needed for finding longest path in the subtree and subtracting its length
 	# Added at 08.10 for testing whether this will improve correspondence between simulation_observed and simulation_expected.
 	sub entrenchment_visitor {
 		my $node = $_[0];
 		my $step = $_[1]->[0];
 		my $subtract_tallest =$_[1]->[2];
 		
	#	my $site_index = $_[1]->[0]; #no site index needed as an argument for this function
	#	my $starting_node = $_[1]->[1];
	#	my %closest_ancestors = %{$_[1]->[2]};
	#	my $depth = $_[2] - $starting_node->get_branch_length ;
	#	if ($node eq $starting_node){return;}
	#	if ($starting_node -> is_terminal){return;}
		if (!$node->is_root){
		my %closest_ancestors = %{$node->get_parent->get_generic("-closest_ancestors")};
		
		if (%closest_ancestors){

		foreach my $site_index(keys %closest_ancestors){ 
			my $anc_node = $closest_ancestors{$site_index};
			my $depth = $static_distance_hash{$anc_node->get_name()}{$node->get_name()};
			#if ($anc_node eq $node) {next;}
			$static_subtree_info{$anc_node->get_name()}{$site_index}{"hash"}{bin($depth,$step)}[1] += $node->get_branch_length;

			my $current_maxdepth = $static_subtree_info{$anc_node->get_name()}{$site_index}{"maxdepth"};
		#	print " current maxdepth $current_maxdepth\n";
			if ($current_maxdepth < $depth){
					$static_subtree_info{$anc_node->get_name()}{$site_index}{"maxdepth"} = $depth;
					if ($subtract_tallest){
							$static_subtree_info{$anc_node->get_name()}{$site_index}{"maxdepth_node"} = $node->get_name();
					}
			}					
		
		}
		
		my @ancestors = keys %closest_ancestors;	
		foreach my $site_index(@ancestors){
			if (!has_no_background_mutation($node, $site_index)){
				delete $closest_ancestors{$site_index};
			}
		}	
		}
		
		foreach my $site_index(keys %{$static_subs_on_node{$node->get_name()}}){
			if ($closest_ancestors{$site_index}){
				my $anc_node = $closest_ancestors{$site_index};
				my $depth = $static_distance_hash{$anc_node->get_name()}{$node->get_name()} - ($node->get_branch_length)/2;
			#	print " ancestor ".$anc_node->get_name(). " node ".$node->get_name()." depth $depth\n";
			#	push $static_subtree_info{$anc_node->get_name()}{$site_index}{"nodes"}, \$node;
				$static_subtree_info{$anc_node->get_name()}{$site_index}{"hash"}{bin($depth,$step)}[0] += 1;
			}
			$closest_ancestors{$site_index} = $node;
		}
		
		$node->set_generic("-closest_ancestors" => \%closest_ancestors);
		}
#!
 	}
   
   
   sub bin {
   	my $depth = $_[0];
   	my $step = $_[1];
   	
   	my $bin = int($depth/$step);
   	if (int($depth/$step) == $depth/$step && $depth != 0){
   		$bin -= 1;
   	}
   	return $bin+1;
   }
    
#sub distances_to_subtree_mutationss {
#	my $tree = parse_tree("C:/Users/weidewind/Documents/CMD/Coevolution/Influenza/Kryazhimsky11/Mock/h1.l.r.newick");
#	my $node = $_[1]; # a node, not a nodename
#   $tree -> prune_tips()
#	$tree->visit_depth_first(
#			-pre_daughter => sub {	my $tnode=shift; 
#									print $tnode -> get_name;
#									print "\n";
#									if ($tnode->get_name eq "alaster3"){
#										$tnode = shift;
#										$tnode = shift;
#									}
#									
#
#		}
#		
#		);
#
#}

sub sieve {
	my @nodes = @{$_[0]};
	my %subs_on_node = %{$_[1]};
	my $site_index = $_[2];
	my $derived = $_[3];
	my $same_derived = $_[4];
	my @sieved;
	
	foreach my $node(@nodes){
		my $t_derived = ${$subs_on_node{$$node->get_name()}}{$site_index}->{"Substitution::derived_allele"};
		#print "Derived: ".$t_derived."\n";
		if (($same_derived & $t_derived eq $derived) || (!$same_derived & $t_derived ne $derived)){
			push @sieved, $node;
		}
	}
	return @sieved;
}

sub dist{
	my $node1 = shift;
	my $node2 = shift;
	my $mrca = $node1->get_mrca($node2);
	my $dist = $node1->get_generic("-mutations_on_path")->Size() +
			   $node2->get_generic("-mutations_on_path")->Size() -
			   $mrca ->get_generic("-mutations_on_path")->Size();

}

sub compute_bitvectors{
	my $tree = $_[0];
	my %mutmap = %{$_[1]};
	my $site_index = $_[2];

	my $root = $tree->get_root();
	$root->set_generic("-mutations_on_path" => Bit::Vector->new(0));
	$root->visit_breadth_first(
				-in => sub{
					my $node=shift;
					if (!$node->is_root()){
						my $pnode = $node->get_parent();
						my $pbit_vect = $pnode -> get_generic("-mutations_on_path");
						my $plength = $pbit_vect->Size();
						my $bit_vect = Bit::Vector->new($plength+1);
						if($plength>0){
							$bit_vect -> Interval_Copy($pbit_vect,0,0,$plength);
						}
						if (exists $mutmap{$node->get_name()}{$site_index}) {
							$bit_vect -> Bit_On($plength)
						};
						$node -> set_generic("-mutations_on_path" => $bit_vect);
						#print ($node->get_name()."\t".$bit_vect->to_Bin()."\n");
					}
				},
				-post => sub{ #all daughters have been processed
					my $node=shift;
				}
			);
}

sub key_for_node_pair{
	#my @pair = sort(@_);
	my @pair = sort (${$_[0]}->get_name,  ${$_[1]}->get_name);	
	return "$pair[0]_$pair[1]";
}

{
	my %pairhash;
	sub divergent_plus{
		print " key ";
		print key_for_node_pair($_[0], $_[1]);
		print " div old ";
		print $pairhash{key_for_node_pair($_[0], $_[1])} -> {"div_count"};
		$pairhash{key_for_node_pair($_[0], $_[1])} -> {"div_count"} +=1; 
		print " div after addition ";
		print $pairhash{key_for_node_pair($_[0], $_[1])} -> {"div_count"};
		print "\n";
	}

	sub convergent_plus{
		$pairhash{key_for_node_pair($_[0], $_[1])} -> {"conv_count"} +=1; 
	}
	
	sub set_distance{
		$pairhash{key_for_node_pair($_[0], $_[1])} -> {"distance"} = $_[2]; 
	}
	
	sub pairhash{
		return %pairhash;
	}
	
	sub value_is_in_array{
		my $value = $_[0];
		my @array = @{$_[1]};
		
		foreach my $v(@array){
			if ($value eq $v){
				return 1;
			}
		}
		return 0;
	}
	

	
}


