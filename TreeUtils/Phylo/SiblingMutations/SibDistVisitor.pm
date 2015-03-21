package TreeUtils::Phylo::SiblingMutations::SibDistVisitor;
#This module provides function for searching sibling mutations
use strict;
use Bio::Phylo::IO;
use Class::Struct;
use Bit::Vector;
use TreeUtils::Phylo::SiblingMutations::Find;


struct SitePairMatrix => {
	site_set => '$', #instance of Bit::Vector. Set of sites where mutations occured at least once.

	#Data matrix representations:
	rh_matrix_val => '$', #for sparse matrices
	ra_matrix_val => '$' #for dence matrices
};

struct SiteDistanceDistr{
	site => '$',
	homo_distr => '@',
	hetero_distr => '@'
};

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
	my $si=0;
			$self->{DATA}->ra_matrix_val=[];
			if(defined $rh_bgr_mutmap){
				foreach my $node($self->{TREE}->get_nodes){
					$name=$node->get_name();
					foreach my $site(@{$rh_bgr_mutmap->{$name}}){
						if(!defined $ra_site2idx->[0]->{$site}){
							$ra_site2idx->[0]->{$site}=$si++;
							push @{$self->{IDX2SITE}->[0]},$site;
						};
					}
				};
				$si=0;
			}
			foreach my $node($self->{TREE}->get_nodes){
				$name=$node->get_name();
				foreach my $site(@{$rh_mutmap->{$name}}){
					if(!defined $ra_site2idx->[$is_intergen]->{$site}){
						$ra_site2idx->[$is_intergen]->{$site}=$si++;
						push @{$self->{IDX2SITE}->[$is_intergen]},$site;
					};
				}
			};
			my $n=@{$self->{IDX2SITE}->[0]};
			$self->{DATA}->site_set1=Bit::Vector->new($n);
			$n=@{$self->{IDX2SITE}->[$is_intergen]};
			$self->{DATA}->site_set2=Bit::Vector->new($n);
			for(my $i=0;$i<@{$self->{IDX2SITE}->[0]};$i++){
				$self->{DATA}->ra_matrix_val->[$i]=[];
				$self->{DATA}->site_set1->Bit_On($i);
				for(my $j=0;$j<@{$self->{IDX2SITE}->[$is_intergen]};$j++){
					$self->{DATA}->ra_matrix_val->[$i]->[$j]=$mtx_init_value;
					$self->{DATA}->site_set2->Bit_On($j);
				}
			}
};
#Methods
sub init_search{
	my $self=shift;
	my $node=shift;
	return if $node->is_root;
	my $dist_mtx=$self->{DATA};
	my $pnode=$node->get_parent;
	$pnode->set_generic("-sites_access" => []);
	my $bvec=$dist_mtx->site_set1->Shadow();
	$bvec->Or($dist_mtx->site_set1,$dist_mtx->site_set2);
	$pnode->get_generic("-sites_access")->[0]=$bvec;
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