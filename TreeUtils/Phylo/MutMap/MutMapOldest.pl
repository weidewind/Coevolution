#!/usr/bin/perl
## This class provides a map of mutations (node_name -> array of sites that mutated on the branch),  
## given a tree and sequences of all nodes (including internal ones)

use strict;
use Bio::Phylo::IO;
use DnaUtilities::compare qw(nsyn_substitutions);

use TreeUtils::Phylo::FigTree;
use TreeUtils::Phylo::MutMap::PhyloDistance;
use Bio::SeqIO;
use Data::Dumper;
use Bit::Vector;
use Statistics::Test::WilcoxonRankSum;
use Try::Tiny;


# creates new object MutMap, which contains tree, sequences of its nodes and the derived map of mutations

sub new{
	my $class=shift;
	my $self = {};
	bless($self,$class);
	
	my $tree = parse_tree($_[0]);
	my %sequences = parse_fasta($_[1]);

	my @ks = keys %fasta;
	my $length = length($sequences{$ks[0]}); # supposing that all sequences have the same length

	$self = {
		tree => $tree,
		sequences =>  %sequences,
		mutmap => mutmap($tree, $sequences),
		seq_length => $length,
	};
	
	return $self;
};


			
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
	return %nodeseqs;
}


## This method gets two files (.newick tree and .fasta with sequences of all its nodes)
## and returns a map of mutations (%subs_on_node: node name -> array of substitution structs; 
## 								   %nodes_with_sub: site index -> array of nodes (links))
sub mutmap {

	my $tree = $_[0];
	my %nodeseqs = %{$_[1]};
	my %subs_on_node;
	my %nodes_with_sub;
	
	my @nodes = $tree -> get_nodes;
	foreach my $node(@nodes){
		if ($node->is_root()) {next;}
		my $name = $node ->get_name();
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



# prints nexus tree, on wich all mutations in the specified site are shown 

sub print_tree_with_mutations{
my $site = shift;
my $tree = $self->{tree};
my @mutmaps = $self->{mutmap};
my %subs_on_node = %{$mutmaps[0]};
my %nodes_with_sub = %{$mutmaps[1]};

my %sites;
my %color;
foreach my $n(@{$nodes_with_sub{$site}}){
	$sites{$$n->get_name()} = $$n->get_name()."_".$site.${$subs_on_node{$$n->get_name()}}{$site}->{"Substitution::derived_allele"};
	$color{$$n->get_name()} = "-16776961";
}

open TREE, ">n1_sites_".$site.".tre";
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




sub logic{
my $tree = $self->{tree};
my $fasta = $self->{sequences};
my @mutmaps = $self->{mutmap};
my %subs_on_node = %{$mutmaps[0]};
my %nodes_with_sub = %{$mutmaps[1]};


my $wilcox_test = Statistics::Test::WilcoxonRankSum->new();
#470
for (my $ind = 1; $ind <566; $ind++){
	compute_bitvectors($tree, \%subs_on_node, $ind);
	my @distr = find_min_distances_naive($ind);

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
 	 try {
 	 	$wilcox_test->load_data(\@{$distr[0]}, \@{$distr[1]});
  		my $prob = $wilcox_test->probability();
  		my $pf = sprintf '%f', $prob; # prints 0.091022
  		print "\n";
  		print $wilcox_test->probability_status();
  		print $wilcox_test->summary();
   		print "\n";
  }
}



}


## 

sub logic_unrestricted{
my $tree = $self->{tree};
my $fasta = $self->{sequences};
my @mutmaps = $self->{mutmap};
my $length  = $self->{seq_length};
my %subs_on_node = %{$mutmaps[0]};
my %nodes_with_sub = %{$mutmaps[1]};

my @general_same;
my @general_diff;

my $wilcox_test = Statistics::Test::WilcoxonRankSum->new();
for (my $ind = 1; $ind <=($length/3); $ind++){

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
	
  	try {
  		$wilcox_test->load_data(\@{$distr[0]}, \@{$distr[1]});
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

# computes minimal_distance_naive for all sites
# returns two arrays: 1) minimal distances between mutations in the same aa,
# 2) min dists between mutations in different aas.

sub logic_general {
	my $tree = parse_tree("C:/Users/weidewind/Documents/CMD/Coevolution/Influenza/Kryazhimsky11/h3.l.r.newick");
my %fasta = parse_fasta("C:/Users/weidewind/Documents/CMD/Coevolution/Influenza/Kryazhimsky11/h3.all.fa");
my @mutmaps = mutmap($tree, \%fasta);
my %subs_on_node = %{$mutmaps[0]};
my %nodes_with_sub = %{$mutmaps[1]};

my @general_same;
my @general_diff;


my $wilcox_test = Statistics::Test::WilcoxonRankSum->new();
#470
for (my $ind = 1; $ind <566; $ind++){
	compute_bitvectors($tree, \%subs_on_node, $ind);
#print ("Nodes with sub: ".scalar @{$nodes_with_sub{$ind}});
my @distr = find_min_distances_naive($tree, \@{$nodes_with_sub{$ind}}, \%subs_on_node, $ind, 1);

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

##control distribution for site i - min distances to any mutation at each site (set size ~ length of the protein*number of branches with mutation in site i)
sub logic_2{
	my $tree = parse_tree("C:/Users/weidewind/Documents/CMD/Coevolution/Influenza/Kryazhimsky11/h1.l.r.newick");
	my %fasta = parse_fasta("C:/Users/weidewind/Documents/CMD/Coevolution/Influenza/Kryazhimsky11/h1.all.fa");
	my @mutmaps = mutmap($tree, \%fasta);
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



## For every node of the given tree sets a hash -min_distances_in_subtree:
## key = index of site (in protein sequence)
## value = minimal distance to any mutation in this site that happened in the subtree down of this node.
## (to find whole-tree minimal distance, you need such arrays for all nodes on the path from this node to the root)


sub compute_min_distances_global{
	my $node = $_[0];
	my %subs_on_node = %{${$self->{mutmap}}[0]};
	
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

# Self-describing. Needs nothing. For each node and for each site computes min distance from that node to 
# closest node in the subtree with mutation in that site.
# Does not care about mutation type. 
# The resulting hash will only contain key "site_index", if there is a mutation in this site somewhere down from that node.
# These hashes are kept in "-min_distances_in_subtree" field of nodes.

sub compute_min_distances_in_subtree{
	my $tree = $self->{tree};
	my %subs_on_node = %{${$self->{mutmap}}[0]};
	
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
					$node->set_generic("-min_distances_in_subtree" => \%min_dists);

					
				}
			);
}


# The only difference from find_min_distances_naive is that this method does not check 
# if there is a substitution on the path between two nodes (if there is a birch between two firs - we still consider this pair of firs)  

sub find_min_distances_unrestricted {
	my $tree = $self->{tree};
	my @nodes = ${$self->{mutmap}}[1]{$_[0]};
	my %subs_on_node = %{${$self->{mutmap}}[0]};
	my $site_index = $_[0];

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

# Given a specific site, for every node with a mutation in this site 
# calculates the minimal distance from this node to another node with a substitution in this site:
# 1) derived aa is the same aa as in the first node 2) derived aa is some other aa
# Does not consider a path between two nodes, if a substitution of another type is present on this path.
# Returns two unsorted arrays of minimal distances.

sub find_min_distances_naive{
	my $tree = $self->{tree};
	my @nodes = ${$self->{mutmap}}[1]{$_[0]};
	my %subs_on_node = %{${$self->{mutmap}}[0]};
	my $site_index = $_[0];

	my @min_distances_same;
	my @min_distances_diff;
	
	for (my $i = 0; $i < scalar @nodes; $i++){
		my $min_dist_same;
		my $min_dist_diff;
		my $bit_vect1 = ${$nodes[$i]}->get_generic("-mutations_on_path")->Clone();
		$bit_vect1->Move_Left(1);
		my $derived1 = ${$subs_on_node{${$nodes[$i]}->get_name()}}{$site_index}->{"Substitution::derived_allele"};
		for (my $j = 0; $j < scalar @nodes; $j++){
			if ($j == $i){ next; }
			my $mrca = ${$nodes[$i]}->get_mrca(${$nodes[$j]})->get_generic("-mutations_on_path")->Size();
			my $bit_vect2 = ${$nodes[$j]}->get_generic("-mutations_on_path")->Clone();
			my $bit_vect1_t = $bit_vect1->Clone();
			$bit_vect2->Move_Left(1);
			$bit_vect2->Move_Right($mrca+1);
			$bit_vect1_t->Move_Right($mrca+1);
			if ($bit_vect2->is_empty() & $bit_vect1_t->is_empty()){
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
	my $tree = $self->{tree};
	my %mutmap = %{$self->{mutmap}};
	my $site_index = $_[0];

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




