#!/usr/bin/perl
use strict;
use Bio::Phylo::IO;

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


