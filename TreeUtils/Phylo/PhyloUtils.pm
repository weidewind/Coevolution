#This package provides utilities for Bio::Phylo::Forest::Tree object 
#to remove zero-length branches and convert them into multifurcations
use strict;
use Bio::Phylo::IO;

sub remove_zero_branches {
		my $tree=shift;
		my @delete;
        $tree->visit_depth_first(
            '-post' => sub {
                my $node = shift;
                my @children = @{ $node->get_children };
                
                #Êthe node is interior, now need to check for each child
                # if it's interior as well
                if ( @children ) {
                    my $has_zero_child;
                    # iterate over children 
                    for my $child ( @children ) {
                        	if ( $child->get_branch_length == 0 ){
                        		$has_zero_child = 1;
                        		last;
                        	}
                    }
                        
                    if ( $has_zero_child ) {
                    	my $parent = $node->get_parent;
                    	for my $child ( @children ) {
                    		my $length = (  $child->get_branch_length || 0 )
                                       + ( $node->get_branch_length || 0 );
                            $child->set_branch_length($length);
                            $child->set_parent($parent);
                    	}
                            $parent->delete($node);
                            # will delete these nodes from the tree array
                            # after the recursion
                            push @delete, $node;						
                        }
                    }				
                }
            
        );
        $tree->delete($_) for @delete;
        return $tree;
}
1;