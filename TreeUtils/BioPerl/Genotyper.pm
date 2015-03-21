#This script determines genotype of query sequences using phylogenetic analysis
#Genotype determines relative to rooted phylogenetic tree of reference sequences
package TreeUtils::BioPerl::Genotyper;
use strict;
use vars qw(@ISA @EXPORT @EXPORT_OK $VERSION);
use Exporter;
$VERSION = 1.00; # Or higher
@ISA = qw(Exporter);
# Symbols to autoexport (:DEFAULT tag)
@EXPORT = qw(genotype_of_query_sequences get_query_atis get_query_genotypes); 
@EXPORT_OK = qw(genotype_of_query_sequences get_query_atis get_query_genotypes);
use Bio::TreeIO;

#Interface
#Calculates CC between annotation and query sequences.
#Use this function to initialize module
sub genotype_of_query_sequences;
#returns hash of sequence ATI (annotation topology indexes). Keys are sequence names.
sub get_query_atis;
#returns sequence genotype string/ Keys  are sequence names
sub get_query_genotypes;



sub annotation_topology_index{
	my ($tree,$ann_mrca,$genotype,$ann_genotype)=@_;
	return undef unless($tree && $ann_mrca && $genotype);
	return 0 unless $ann_genotype;
	return 0 if $ann_genotype==$ann_mrca;
	return 1 if $ann_genotype==$genotype;
	my $d1=$tree->distance(-nodes => [($ann_mrca,$ann_genotype)]);
	my $d2=$tree->distance(-nodes => [($ann_mrca,$genotype)]);
	return $d2/$d1;
};	

sub topology_index{
	my ($tree,$outgroup_id,$genotype)=@_;
	my $root_node=$tree->get_root_node;
	return undef unless($tree && $genotype);
	return 0 if $genotype==$root_node;
	my ($Node)=$tree->find_node($outgroup_id);
	my $node;
	foreach my $inode($root_node->each_Descendent){
		next if $inode==$Node;
		$node=$inode;
		last;
	};
	my $d1=$tree->subtree_length($node);
	my $d2=$tree->subtree_length($genotype);
	my $d3=$tree->distance(-nodes => [($node,$genotype)]);
	return ($d1-$d2-$d3)/$d1;
};	

#genotypes of query sequences
my %query_genotypes=();
#query annotation topology indexes
my %query_atis=();


#ARGS:
#rt_file [in] - A file with a reference tree
#ann_file [in] - A FASTA list of tips in a reference tree used for calculation topological index (Not need to align. Only sequence names used!)
#qt_file [in] - A query tree - a tree with tips corresponded to reference and query sequences
#queries [in] - A FASTA file with query sequences (Not need to align. Only sequence names used!)
#outgroup_id [in] - A sequence name used as outgroup
sub genotype_of_query_sequences{
	my ($rt_file,$ann_file,$qt_file,$query_file,$outgroup_id)=@_;
	my $ref_treeio = Bio::TreeIO->new(-format => 'newick', -file => $rt_file);
	my $ref_tree=$ref_treeio->next_tree;
	if(!$ref_tree){
		warn "\nWrong referendce tree file: $rt_file";
		return undef;
	};
	if(!$outgroup_id){
		warn "\nWrong outgroup: $outgroup_id";
		return undef;
	};
	my ($root_node)=$ref_tree->find_node($outgroup_id);
	if(!$root_node){
		warn "\nA sequence with id: $outgroup_id doesn't exist in the reference tree!";
		return undef;
	};
	if(!$root_node->ancestor){
		warn "\nA leaf node: $outgroup_id doesn't have an ancestor in the reference tree!";
		return undef;
	};
	$ref_tree->reroot_at_midpoint($root_node);
	my $qt_treeio = Bio::TreeIO->new(-format => 'newick', -file => $qt_file);
	my $qt_tree;
	if(!($qt_tree = $qt_treeio->next_tree)){
		warn "\nWrong query tree file: $qt_file" ;
		return undef;
	};
	($root_node)=$qt_tree->find_node($outgroup_id);
	if(!$root_node){
		warn "\nA sequence with id: $outgroup_id doesn't exist in the query tree!";
		return undef;
	};
	if(!$root_node->ancestor){
		warn "\nA leaf node: $outgroup_id doesn't have an ancestor in the query tree!";
		return undef;
	};
	$qt_tree->reroot_at_midpoint($root_node);
	#read set of queries
	open FILE, $query_file or die "\nUnable to open file with queries: $query_file!";
	my %query_ids;
	my %query_weights;
	my $W=0;
	while(<FILE>){
		if(/\s*>\s*(\S+)/){
			my $qid=$1;
			$query_ids{$qid}=1;
			if(/\s*>\s*\S+\s+([0\.\d]+)/){
				$query_weights{$qid}=$1;
				$W+=$query_weights{$qid};
			};
		};
	};
	foreach my $qid(keys %query_weights){
		$query_weights{$qid}/=$W;
	};
	close FILE;
	my $IsQueryWeghted=1;
	if(scalar(keys %query_ids)!=scalar(keys %query_weights)){
		%query_weights=();
		warn "\nNot all sequences in the file: $query_file\n\thave weights!";
		$IsQueryWeghted=0;
	};
	#read set of queries
	open FILE, $ann_file or die "\nUnable to open file with queries: $ann_file!";
	my %ann_ids;
	while(<FILE>){
		if(/\s*>\s*(\S+)/){
			$ann_ids{$1}=1;
		};
	};
	close FILE;
	
	#loocking for a commom ancestor of annotating sequences on the reference tree
	my $ann_mrca;
	my @ann_nodes;
	foreach my $aid(keys(%ann_ids)){
		my ($anode)=$ref_tree->find_node($aid);
		if(!$anode){
			warn "\nThe annotating sequence: $aid doesn.t exist in the reference tree!";
			return undef;
		};
		push @ann_nodes,$anode;
	};
	$ann_mrca=$ref_tree->get_lca(@ann_nodes);
	
	#Genotype of a query sequence is a MRCA of all it's most close references
	my $CC=0;
	my $N=0;
	foreach my $qid(keys %query_ids){
		my ($qnode)=$qt_tree->find_node($qid);
		my $anc=$qnode->ancestor;
		next unless $anc;
		my @ref_family;
		my @qt_family;   
		my $inode=$qnode;
		my $str=0;
		my $n=1;
		#look for a most close annotation reference
		while($str==0||$n<3){
			my @qtnodes;
			my $rf_vol=@ref_family;
			#get all leaf under $anc node
			foreach my $node1($anc->each_Descendent){
				next if $node1==$inode; #do not get down into already processed subtree 
				if(!$node1->is_Leaf){
					foreach my $node2($node1->get_all_Descendents){
						push @qtnodes, $node2 if($node2->is_Leaf);
					};
				}else{
					push @qtnodes, $node1;
				};
			};
			foreach my $qnode(@qtnodes){
				my $id=$qnode->id;
				next if $query_ids{$id};
				$str++ if $ann_ids{$id}; #unnotation found!
				my ($rnode)=$ref_tree->find_node($id);
				if(!$rnode){
					warn "\nA reference tree doesn't have reference tip with id=$id";
					return undef;
				}; 
				push @ref_family, $rnode;
			};
			$inode=$anc;
			$anc=$anc->ancestor;
			$n++ if($rf_vol<@ref_family);
			last unless $anc;
		};
		if($str==0){
			warn "\nUnable to find annotation for qid=$qid!";
			return undef;
		};
		
		#determing genotype (mrca) of obtained reference homologs using the reference tree
		if(@ref_family>1){
			$anc=$ref_tree->get_lca(@ref_family);
		}else{
			$anc=$ref_family[0]->ancestor;
		};
		next if(!$anc);
		@ref_family=();
		@qt_family=();
		#add all references correcponded to $anc into a node family
		foreach my $rnode($anc->get_all_Descendents){
			if($rnode->is_Leaf){
				my $id=$rnode->id;
				my ($qnode)=$qt_tree->find_node($id);  
				if(!$qnode){
					warn "\nA query tree doesn't have reference tip with id=$id";
					return undef;
				};
				push @qt_family, $qnode;
			};
		};
		if(@qt_family>1){
			$anc=$qt_tree->get_lca(@qt_family);
		}else{
			$anc=$qt_family[0]->ancestor;
		};
		foreach my $qnode($anc->get_all_Descendents){
			if($qnode->is_Leaf){
				my $id=$qnode->id;
				next if $query_ids{$id};
				my ($rnode)=$ref_tree->find_node($id);  
				if(!$rnode){
					warn "\nA reference tree doesn't have reference tip with id=$id";
					return undef;
				};
				push @ref_family, $rnode;
			};
		};
		#determs a genotype of query sequence
		if(@ref_family>1){
			$anc=$ref_tree->get_lca(@ref_family);
		}else{
			$anc=$ref_family[0]->ancestor;
		};
		my @nodes;
		foreach my $rnode(@ref_family){
			my $id=$rnode->id;
			next unless $ann_ids{$id};
			push @nodes, $rnode;
		};
		$inode=undef;
		if(@nodes==1){
			$inode=$nodes[0]->ancestor;
		}elsif(@nodes>1){
			$inode=$ref_tree->get_lca(@nodes);
		};
		my $id=$ref_family[0]->id();
		$str=$id;
		#print "\n$qid\t($id";
		for(my $i=1;$i<@ref_family;$i++){
			$id=$ref_family[$i]->id;
			$str.=" $id";
			#print ", $id";
		};
		#print ")";
		$query_genotypes{$qid}=$str;
		my $ati=annotation_topology_index($ref_tree,$ann_mrca,$anc,$inode);
		$query_atis{$qid}=$ati;
		#my $ti=topology_index($ref_tree,$outgroup_id,$anc);
		#print "\n$qid\t$ati";
		$N++;
		if($IsQueryWeghted){
			$CC+=$query_weights{$qid}*$ati;
		}else{
			$CC+=$ati;
		};
	};
	if(!$IsQueryWeghted){
		$CC/=$N;
	};
	return $CC;
};

sub get_query_genotypes{
	return %query_genotypes;
};

sub get_query_atis{
	return %query_atis;
};

1;