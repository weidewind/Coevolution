package TreeUtils::Phylo::SiblingMutations::BasicVisitor;
#This module provides function for searching sibling mutations
use strict;

#interface declaration
#constructor
sub new{
	my $class=shift;
	my $self={};
	bless($self,$class);
	$self->_init(@_);
	return $self;
};

sub _init{
	my $self=shift;
	if(@_){
		my %args = @_;
		while(($k,$v)=each %args){
			$self->{$k}=$v;
		}
	};
};
#Methods
sub init_search{
	my $self=shift;
	my $node=shift;
	return if $node->is_root;
	my $str=$node->get_name;
	print "Traversing prallel lineages for the node: $str\nMSA\tbranch\tdistance";
	my $pnode=$node->get_parent();
	$pnode->set_score(0);
};

sub clean_node{
	my $self=shift;
	my $node=shift;
	$node->set_generic();
	$node->set_score(0);
};

#Can change visitor's internal data. Does not change nodes!
sub update_visitor{
	my $self=shift;
	my ($src_node,$mca_node,$node)=@_;
	my $l=$pnode->get_score()+$node->get_branch_length()/2;
	$l-=2*$mca_node->get_score()+$src_node->get_branch_length()/2;
	$l*-1 if $l<0;
	my $str=$mca_node->get_name();
	print "\n$str";
	$str=$node->get_name();
	print "\t$str\t$l";
	return 1;
};

sub update_node{
	my $self=shift;
	my ($node_source,$node_branch,$node_target)=@_;
	my $l=$node_source->get_score()+$node_branch->get_branch_length();
	$node_target->set_score($l);
	return 1;
};

1;