package TreeUtils::Phylo::SiblingMutations::ExpStatistics;
#This module calculates exponential distance statistics for sibling mutations
@ISA = qw(TreeUtils::Phylo::SiblingMutations::BasicVisitor);

use strict;
use Bio::Phylo::IO;
use Class::Struct;
use Bit::Vector;
use TreeUtils::Phylo::SiblingMutations::Find;

#matrix of phylogenetic distances between sites
my $mtx_init_value=0;
my $tau=50;

struct SitePairMatrix => {
	site_set1 => '$', #instance of Bit::Vector
	site_set2 => '$', #instance of Bit::Vector
	gene_index1 => '$',
	gene_index2 => '$', #gene_index2 == gene_index1 in the case of both sites within the same gene
	#Data matrix representations:
	rh_matrix_val => '$', #for sparse matrices
	ra_matrix_val => '$' #for dence matrices
};

struct SitePairExpStatistics{
	bg_site => '$',
	tg_site => '$',
	value => '$'
};

#index access
sub _site_pair{
	my $self=shift;
	if(@_){
		my ($i,$j,$v)=@_;
		my $mtx=$self->{DATA};
		my $tmp;
		if(defined $mtx->rh_matrix_val){
			$tmp=$mtx->rh_matrix_val->{"$i,$j"};
			if(defined $tmp){
				$mtx->rh_matrix_val->{"$i,$j"}=$v if defined $v;
			};
		}elsif(defined $mtx->ra_matrix_val){
			$tmp=$mtx->ra_matrix_val->[$i]->[$j];
			$mtx->ra_matrix_val->[$i]->[$j]=$v if defined $v;
		};
		return $tmp;
	};
	return undef;
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
		my ($in_tree,$rh_mutmap,$rh_bgr_mutmap,$ra_site_pairs);
		my %args = @_;
		$self->{TAU}=$tau;
		while(($k,$v)=each %args){
			if($k eq "-tree"){
				$in_tree=$v;
			}elsif($k eq "-mutation_map"){
				$rh_mutmap=$v;
			}elsif($k eq "-bgr_mutation_map"){
				$rh_bgr_mutmap=$v;
			}elsif($k eq "-site_pairs"){
				$ra_site_pairs=$v;
			}elsif($k eq "-tau"){
				$$self->{TAU}=$v;
			}else{
				die "\nUnknown parameter: $k!";
			};
			$self->{$k}=$v;
		};
		die "\nInput tree was not specified!" unless defined $in_tree;
		die "\nNot defined mutations on the input tree branches!" unless defined $rh_mutmap;
		$self->{TREE}=$in_tree->clone();
		$self->{TREE}->set_generic();
		my $is_intergen=0; 
		$is_intergen=1 if(defined $rh_bgr_mutmap); #f: directional mutation in one gene, $is_intergen=0
		my $ra_site2idx=[];
		$ra_site2idx->[0]={}; #f
		$ra_site2idx->[1]={} if $is_intergen;
		$self->{IDX2SITE}=[];
		$self->{IDX2SITE}->[0]=[]; #f
		$self->{IDX2SITE}->[1]=[] if $is_intergen;
		$self->{DATA}=SitePairMatrix->new();
		$self->{DATA}->gene_index1=0;
		if(defined $rh_bgr_mutmap){
			$self->{DATA}->gene_index2=1;
		}else{
			$self->{DATA}->gene_index2=0;
		};
		if(defined $ra_site_pairs){
			$self->{DATA}->rh_matrix_val={};
			my @si=(0,0); #site counters for two genes
			my %sites1;
			my %sites2;
			foreach my $spi(@{$ra_site_pairs}){ # индексируем все сайты в выданных нам парах?
				if(!defined $ra_site2idx->[0]->{$spi->bg_site}){
					$ra_site2idx->[0]->{$spi->bg_site}=$si[0]++;
					push @{$self->{IDX2SITE}->[0]},$spi->bg_site;
				};
				if(!defined $ra_site2idx->[$is_intergen]->{$spi->tg_site}){
					$ra_site2idx->[$is_intergen]->{$spi->tg_site}=$si[$is_intergen]++;
					push @{$self->{IDX2SITE}->[$is_intergen]},$spi->tg_site;
				};
				my ($s1,$s2)=($ra_site2idx->[0]->{$spi->bg_site},$ra_site2idx->[$is_intergen]->{$spi->tg_site});
				
			#в хэше sites1  ключи - индексы всех возможных сайтов бэкграундных мутаций  
				$sites1{$s1}=1;
				$sites2{$s2}=1;
				$self->{DATA}->rh_matrix_val->{"$s1,$s2"}=$mtx_init_value;
			};
			my $n=@{$self->{IDX2SITE}->[0]};
			my @tmp=keys %sites1;
			$self->{DATA}->site_set1=Bit::Vector->new($n);
			foreach my $i(@tmp){
				$self->{DATA}->site_set1->Bit_On($i); # вот это что? мне казалось, что индексы в бит-векторе должны соответствовать расстоянию от ветки до корня?
													  #	ап, это вообще не те бит-векторы) Просто сет сайтов, в которых происходили мутации. Наверное.
			};
			$n=@{$self->{IDX2SITE}->[1]} if $is_intergen;
			@tmp=keys %sites2;
			$self->{DATA}->site_set2=Bit::Vector->new($n);
			foreach my $i(@tmp){
				$self->{DATA}->site_set2->Bit_On($i);
			};
		}else{#no restriction for site pairs, use dence matrix data representation
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
		foreach my $node($self->{TREE}->get_nodes){
			$node->set_generic("-sites_on_branch" => []);
			$node->get_generic("-sites_on_branch")->[0]=[];
			$node->get_generic("-sites_on_branch")->[1]=[] if $is_intergen;
			$name=$node->get_name();
			if(defined $rh_bgr_mutmap){
				foreach my $site(@{$rh_bgr_mutmap->{$name}}){
					my $sidx=$ra_site2idx->[0]->{$site};
					push @{$node->get_generic("-sites_on_branch")->[0]},$sidx if(defined $sidx);
				};
			};
			foreach my $site(@{$rh_mutmap->{$name}}){
				my $sidx=$ra_site2idx->[$is_intergen]->{$site};
				push @{$node->get_generic("-sites_on_branch")->[$is_intergen]},$sidx if(defined $sidx);
			};
		}
	};
	$self->{SITE2IDX}=$ra_site2idx;
	$self->_calculate();
};

sub gen_matrix_iterator{
	my $self=shift;
	my ($s1_idx,$s2_idx);
	my ($s1_last,$s2_last);
	my $is_intergen=@{$self->{IDX2SITE}}-1;
	my @mtx_keys;
	if(defined $self->{DATA}->rh_matrix_val){
		@mtx_keys=keys %{$self->{DATA}->rh_matrix_val};
		$s1_idx=0;
		$s1_last=@mtx_keys;
	}elsif(defined $self->{DATA}->ra_matrix_val){
		$s1_idx=0;
		$s1_last=@{$self->{IDX2SITE}->[0]};
		$s2_idx=0;
		$s2_last=@{$self->{IDX2SITE}->[$is_intergen]};
	};
	return sub {
		my $v=SitePairExpStatistics->new();
		if(defined $self->{DATA}->rh_matrix_val){
			return undef if $s1_idx == $s1_last;
			my @idx=split ',', $mtx_keys[$s1_idx];
			$v->bg_site($self->{DATA}->{IDX2SITE}->[0]->[$idx[0]]);
			$v->tg_site($self->{DATA}->{IDX2SITE}->[$is_intergen]->[$idx[1]]);
			$v->value($self->{DATA}->rh_matrix_val->{$mtx_keys[$s1_idx]});
			$s1_idx++;
			return $v; 
		}elsif(defined $self->{DATA}->ra_matrix_val){
			return undef if $s1_idx == $s1_last || $s2_idx == $s2_last;
			$v->bg_site($self->{DATA}->{IDX2SITE}->[0]->[$s1_idx]);
			$v->tg_site($self->{DATA}->{IDX2SITE}->[$is_intergen]->[$s2_idx]);
			$v->value($self->{DATA}->ra_matrix_val->[$s1_idx]->[$s2_idx]);
			$s2_idx++;
			if($s2_idx==$s2_last){
				$s1_idx++;
				$s2_idx=0;
			};
			return $v;
		};
		return undef;
	};
};
#Private

sub _calculate{
	my $self=shift;
	foreach my $node(@{$self->{TREE}->get_nodes}){
		next if $node->is_root;
		findSiblings($node,$self);
	}
};

sub sites_on_branch2bitvect{
	my ($node,$gene_idx,$bvect_shadow)=@_;
	foreach my $site_idx(@{$node->get_generic("-sites_on_branch")->[$gene_idx]}){
		$bvect_shadow->Bit_On($idx);
	}
};

sub update_node{
	my $self=shift;
	my ($node_source,$node_branch,$node_target)=@_;
	$node_target->set_generic("-sites_access" => []);
	my $s_ok=1;
	for(my $idx=0;$idx<@{$node_source->get_generic("-sites_access")};$idx++){
		my $tgt_bvec=$node_source->get_generic("-sites_access")}->[$idx]->Clone();
		foreach my $site_idx(@{$node_branch->get_generic("-sites_on_branch")->[$idx]}){
			$tgt_bvec->Bit_Off($site_idx);
		};
		$node_target->get_generic("-sites_access")->[$idx]=$tgt_bvec;
		$s_ok=0 if $tgt_bvec->is_empty();
	};
	my $l=$node_source->get_score()+$node_branch->get_branch_length();
	$node_target->set_score($l);
	return $s_ok;
};

sub update_visitor{
	my $self=shift;
	my $site_pair_mtx=$self->{DATA};
	my ($src_node,$mca_node,$node)=@_;
	my $pnode=$node->get_parent;
	my $sites1_bvec=$pnode->get_generic("-sites_access")->[$site_pair_mtx->gene_index1];
	my $sites2_bvec=$pnode->get_generic("-sites_access")->[$site_pair_mtx->gene_index2];
	return 0 if $sites1_bvec->is_empty() || $sites2_bvec->is_empty();
	my $both_sites1;
	my $both_sites2;
	{
		#get bit vectors on source and target branches
		my $src_sites1=$sites1_bvec->Shadow();
		sites_on_branch2bitvect($src_node,$site_pair_mtx->gene_index1,$src_sites1);
		my $tgt_sites1=$sites1_bvec->Shadow();
		sites_on_branch2bitvect($node,$site_pair_mtx->gene_index1,$tgt_sites1);
		$both_sites1=$src_sites1->Shadow();
		$both_sites1->And($src_sites1,$tgt_sites1);
		if($site_pair_mtx->gene_index1 == $site_pair_mtx->gene_index2){
			$both_sites2=$both_sites1;
		}else{
			my $src_sites2=$sites2_bvec->Shadow();
			sites_on_branch2bitvect($src_node,$site_pair_mtx->gene_index2,$src_sites2);
			my $tgt_sites2=$sites2_bvec->Shadow();
			sites_on_branch2bitvect($node,$site_pair_mtx->gene_index2,$tgt_sites2);
			$both_sites2=$src_sites2->Shadow();
			$both_sites2->And($src_sites2,$tgt_sites2);
		}
	};
	my $l=$pnode->get_score()+$node->get_branch_length()/2;
	$l-=2*$mca_node->get_score()+$src_node->get_branch_length()/2;
	$l*-1 if $l<0;
	my $n=0;
	# обновляет данные для каждой пары сайтов, мутации в которых присут на этих ветках.
	foreach my $s2_idx(@{$node->get_generic("-sites_on_branch")->[$site_pair_mtx->gene_index2]}){
		next unless $sites2_bvec->bit_test($s2_idx);
		my $p=1.0;
		$p=0.5 if $both_sites2->bit_test($s2_idx);
		foreach my $s1_idx(@{$src_node->get_generic("-sites_on_branch")->[$site_pair_mtx->gene_index1]}){
			next unless $sites1_bvec->bit_test($s1_idx);
			$p*=0.5 if $both_sites1->bit_test($s1_idx);
			my $v=$self->_site_pair($s1_idx,$s2_idx);
			if(defined $v){
				$v+=$p*exp(-$l/$tau);
				$self->_site_pair($s1_idx,$s2_idx,$v);
				$n++;
			};
		}
	}
	return $n;
};

sub init_search{
	my $self=shift;
	my $node=shift;
	return if $node->is_root;
	my $dist_mtx=$self->{DATA};
	my $pnode=$node->get_parent;
	$pnode->set_generic("-sites_access" => []);
	if($dist_mtx->gene_index1==$dist_mtx->gene_index2){
		my $bvec=$dist_mtx->site_set1->Shadow();
		$bvec->Or($dist_mtx->site_set1,$dist_mtx->site_set2);
		$pnode->get_generic("-sites_access")->[0]=$bvec;
	}else{
		$pnode->get_generic("-sites_access")->[0]=$dist_mtx->site_set1->Clone();
		$pnode->get_generic("-sites_access")->[1]=$dist_mtx->site_set2->Clone();
	};
	$pnode->set_score(0);
};

sub cleanNode{
	my $node=shift;
	$node->set_generic();
	$node->set_score(0);
};

1;