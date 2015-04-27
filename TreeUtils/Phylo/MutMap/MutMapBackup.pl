#!/usr/bin/perl
## This class provides a map of mutations (node_name -> array of sites that mutated on the branch),  
## given a tree and sequences of all nodes (including internal ones)

use strict;
use Bio::Phylo::IO;
use DnaUtilities::compare qw(nsyn_substitutions syn_substitutions);

use TreeUtils::Phylo::FigTree;
use Bio::SeqIO;
use Data::Dumper;
use Bit::Vector;
use Statistics::Test::WilcoxonRankSum;
use Try::Tiny;

use List::Util qw/shuffle/; 
use Statistics::Basic qw(:all);

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

sub synmutmap {

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


sub load_data{
	
}

#print_tree_with_mutations(6);
#print_tree_with_mutations(465);
# prints nexus tree, on wich all mutations in the specified site are shown 

sub print_tree_with_mutations{
my $site = shift;
my $tree = parse_tree("C:/Users/weidewind/Documents/CMD/Coevolution/Influenza/Kryazhimsky11/n1.l.r.newick");
my %fasta = parse_fasta("C:/Users/weidewind/Documents/CMD/Coevolution/Influenza/Kryazhimsky11/n1.all.fa");
my @mutmaps = synmutmap($tree, \%fasta);
my %subs_on_node = %{$mutmaps[0]};
my %nodes_with_sub = %{$mutmaps[1]};

my %sites;
my %color;
foreach my $n(@{$nodes_with_sub{$site}}){
	$sites{$$n->get_name()} = $$n->get_name()."_".$site.${$subs_on_node{$$n->get_name()}}{$site}->{"Substitution::derived_allele"};
	$color{$$n->get_name()} = "-16776961";
}

open TREE, ">n1_syn_sites_".$site.".tre";
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

#onemtest(248);
#onemtest(214);
#onemtest(136);

#logic_general();
#logic_unrestricted();
#logic();

sub logic{
my $tree = parse_tree("C:/Users/weidewind/Documents/CMD/Coevolution/Influenza/Kryazhimsky11/n1.l.r.newick");
my %fasta = parse_fasta("C:/Users/weidewind/Documents/CMD/Coevolution/Influenza/Kryazhimsky11/n1.all.fa");
my @mutmaps = synmutmap($tree, \%fasta);
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
	my $tree = parse_tree("C:/Users/weidewind/Documents/CMD/Coevolution/Influenza/Kryazhimsky11/n2.l.r.newick");
my %fasta = parse_fasta("C:/Users/weidewind/Documents/CMD/Coevolution/Influenza/Kryazhimsky11/n2.all.fa");
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

subtreetest();

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




