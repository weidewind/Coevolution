 use Bio::Phylo::IO;
 use Bio::Phylo::Forest::NodeRole;
 use Bio::Phylo::Forest::TreeRole;
 

sub parseTree {
					my $tree_file = $_[0];
					open TREE, "< $tree_file" or die "Невозможно открыть файл ".$tree_file."\n";
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

#Algorithm performs the "Down And Up" search for all parallel lineages for a specified tree node
#Arguments:
#	$src_node - A tree node. Expects Bio::Phylo interface
#	$visitor - An object that processes algorithm events

sub findSiblings{
	my $src_node=shift;
	my $visitor=shift;
	
	return if $src_node->is_root;
	my $is_ok=0;
	my $node=$src_node;
	my $pnode=$node->get_parent;
	$visitor->init_search($node);
	while(!$node->is_root){
		#go down on each parallel lineage
		foreach my $chnode(@{$pnode->get_children}){
			next if $chnode == $node;
			$chnode->visit_breadth_first(
				-in => sub{
					my $node=shift;
					$is_ok+=$visitor->update_visitor($src_node,$pnode,$node);
					$visitor->update_node($node->get_parent,$node,$node);
				},
				-post => sub{ #all daughters have been processed
					my $node=shift;
					$visitor->clean_node($node);
				}
			);
		};
		#go up
		$node=$pnode;
		$pnode=$node->get_parent;
		if(defined $pnode){
			last unless $visitor->update_node($node,$node,$pnode);
		};
		$visitor->clean_node($node);
	}
	return $is_ok;
};

my $tree = parseTree("C:/Users/weidewind/Documents/CMD/BCV_project/OurSamples/Endocarditis/1354_1423578065169/Output/3655-785R/cluster_1_0.571.mltree1.with_taxonomy");
print $tree -> get_root()