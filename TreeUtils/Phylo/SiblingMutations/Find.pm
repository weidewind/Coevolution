package TreeUtils::Phylo::SiblingMutations::Find;
#This module provides function for searching sibling mutations
use vars qw(@ISA @EXPORT @EXPORT_OK $VERSION);
use Exporter;
$VERSION = 1.00; # Or higher
@ISA = qw(Exporter);
@EXPORT = qw(findSiblings); # Symbols to autoexport (:DEFAULT tag)
@EXPORT_OK = qw(findSiblings);
use strict;
use Bio::Phylo::IO;

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

1;