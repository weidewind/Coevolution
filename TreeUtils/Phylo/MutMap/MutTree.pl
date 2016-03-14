#!/usr/bin/perl
## This class provides a tree, map of mutations (node_name -> array of sites that mutated on the branch),  
## given a tree and sequences of all nodes (including internal ones)

use strict;
use Bio::Phylo::IO;
use DnaUtilities::compare qw(nsyn_substitutions);

use TreeUtils::Phylo::FigTree;
use TreeUtils::Phylo::MutMap::PhyloDistance;