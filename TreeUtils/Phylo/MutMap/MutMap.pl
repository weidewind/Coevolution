#!/usr/bin/perl
## This class provides a map of mutations (node_name -> array of sites that mutated on the branch),  
## given a tree and sequences of all nodes (including internal ones)

use strict;
use Bio::Phylo::IO;
use DnaUtilities::compare qw(nsyn_substitutions);

use TreeUtils::Phylo::FigTree;
use Bio::SeqIO;
use Data::Dumper;
use Bit::Vector;
use Statistics::Test::WilcoxonRankSum;
use Try::Tiny;

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
	return %nodeseqs;
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


sub onemtest{
	my $site = shift;
	my $tree = parse_tree("C:/Users/weidewind/Documents/CMD/Coevolution/Influenza/Kryazhimsky11/n1.l.r.newick");
my %fasta = parse_fasta("C:/Users/weidewind/Documents/CMD/Coevolution/Influenza/Kryazhimsky11/n1.all.fa");
my @mutmaps = mutmap($tree, \%fasta);
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

onemtest(248);
onemtest(214);
onemtest(136);


sub logic{
my $tree = parse_tree("C:/Users/weidewind/Documents/CMD/Coevolution/Influenza/Kryazhimsky11/h1.l.r.newick");
my %fasta = parse_fasta("C:/Users/weidewind/Documents/CMD/Coevolution/Influenza/Kryazhimsky11/h1.all.fa");
my @mutmaps = mutmap($tree, \%fasta);
my %subs_on_node = %{$mutmaps[0]};
my %nodes_with_sub = %{$mutmaps[1]};


#my @sieved = sieve(\@{$nodes_with_sub{10}},\%subs_on_node, 10, "Y", 1);
#print "Sieved: ";
#foreach my $s(@sieved){
#	print $$s->get_name()."\t";
#}
my $wilcox_test = Statistics::Test::WilcoxonRankSum->new();
#470
for (my $ind = 1; $ind <566; $ind++){
	compute_bitvectors($tree, \%subs_on_node, $ind);
#print ("Nodes with sub: ".scalar @{$nodes_with_sub{$ind}});
my @distr = find_min_distances_naive($tree, \@{$nodes_with_sub{$ind}}, \%subs_on_node, $ind, 1);

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

sub subtreetest{
	my $tree = parse_tree("C:/Users/weidewind/Documents/CMD/Coevolution/Influenza/Kryazhimsky11/h1.l.r.newick");
my %fasta = parse_fasta("C:/Users/weidewind/Documents/CMD/Coevolution/Influenza/Kryazhimsky11/h1.all.fa");
my @mutmaps = mutmap($tree, \%fasta);
my %subs_on_node = %{$mutmaps[0]};
my %nodes_with_sub = %{$mutmaps[1]};

 compute_min_distances_in_subtree($tree, \%subs_on_node);
 my $root = $tree->get_root();
	$root->visit_depth_first(
				-post => sub{ #all daughters have been processed
					my $node=shift;
					my %min_dists = compute_min_distances_global($node, \%subs_on_node);
					print $node->get_name."\t";
					for my $site(keys %min_dists){
						print ($site."_".$min_dists{$site}."\t");
					}
					print "\n";
				}
	);
}

#subtreetest();

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




