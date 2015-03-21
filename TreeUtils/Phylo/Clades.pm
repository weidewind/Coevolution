#This package provides utilities for clade identification and statistic calculations on the Bio::Phylo::Forest::Tree object
#The functions in the package assume that a tree object has names for internal nodes
use strict;
use Bio::TreeIO;
use Bio::Phylo::IO;
use Class::Struct;

my @reassortmentMapMode=("default","force_mrca","recent_mrca","mrca_only");
my $rmap_mode_idx=0;

sub setReassortmentMapMode{
	my $mode=shift;
	return "$rmap_mode_idx:$reassortmentMapMode[$rmap_mode_idx]" unless defined $mode;
	my $tmp=$rmap_mode_idx;
	$mode=lc $mode;
	if($mode eq "force_mrca"){
		$rmap_mode_idx=1;
	}elsif($mode eq "recent_mrca"){
		$rmap_mode_idx=2;
	}elsif($mode eq "mrca_only"){
		$rmap_mode_idx=3;
	}else{
		$rmap_mode_idx=0;
	};
	return "$rmap_mode_idx:$reassortmentMapMode[$rmap_mode_idx]";
};
	

#This function finds a set of branches corresponded to the predefined taxon subset @{$ra_tax_names}
sub find_clades{
	my ($tree,$ra_tax_names)=@_;
	my $outgrNode=get_outgroup($tree);
	my $t=1;
	my @clades;
	my %tax_names;
	foreach my $taxn(@{$ra_tax_names}){
		if($t && ($taxn eq $outgrNode->get_name)){
			$t=0;
		};
		$tax_names{$taxn}=1;
	};
	#return () if @{$tree->get_terminals()}/3<@{$ra_tax_names};
	my @taxa;
	foreach my $node(@{$tree->get_terminals()}){
		if(defined($tax_names{$node->get_name})){
			push @taxa,$node if $t;
		}else{
			push @taxa,$node unless $t;
		}
	};
	my @tmp=@taxa;
	my $mrca=$tree->get_mrca(\@tmp);
	@clades=split_clade($mrca,\@taxa);
	if(@clades){
		if($rmap_mode_idx==1){
			push @clades, $mrca;
		}elsif($rmap_mode_idx==2){
			my $dist=0;
			for(my $i=0;$i<@clades-1;$i++){
				for(my $j=$i+1;$j<@clades;$j++){
					my $d=$clades[$i]->calc_patristic_distance($clades[$j]);
					$dist=$d if($d>$dist);
				}
			}
			my $avg_branch_length=0;
			foreach my $node(@{$tree->get_internals()}){
				$avg_branch_length+=$node->get_branch_length if $node->get_branch_length>0;
			};
			$avg_branch_length/=@{$tree->get_internals()};
			push @clades, $mrca if($dist<=4*$avg_branch_length);
		}elsif($rmap_mode_idx==3){
			@clades=();
			push @clades, $mrca;
		}
	}else{
		push @clades, $mrca;
	};
	return @clades;
};

struct Reassortment =>{
	node => '$',
	info => '@',
	nobs => '$',
	conf => '$'
};
struct ReassortmentInfo =>{
	id => '$',
	dataset =>'$',
	conf => '$',
	cardinality => '$'
};
sub new_reassortment{
	my $node=shift;
	my $rinf=Reassortment->new();
	$rinf->node($node);
	$rinf->nobs(1);
	$rinf->conf(0);
	return $rinf;
};

sub read_reassortment_info{
	my $ra_report_files=shift;
	my $tree=shift;
	my $min_rconf=0;
	$min_rconf=$_[0] if defined $_[0];
	my $min_robs=1;
	$min_robs=$_[1] if defined $_[1];
	my $dataset_name=$_[2];
	my $ra_rmtaxa=$_[3];
	my %removed_taxa;
	if(defined($ra_rmtaxa)&&(@{$ra_rmtaxa}>0)){
		foreach my $str(@{$ra_rmtaxa}){
			$removed_taxa{$str}=1;
		}
	}
	my %reassortments;
	foreach my $report_fn(@{$ra_report_files}){
		open REPORT, "<$report_fn" or die "\nUnable to open input file: $report_fn!";
		my $seg_name=$report_fn;
		if($report_fn=~/([^\/]+)_report$/){
			$seg_name=$1;
			if(defined $dataset_name){
				$seg_name=~s/$dataset_name//;
				$seg_name=~s/^-+//;
				$seg_name=~s/-+$//;
			};
		};
		while(<REPORT>){
			if(/Candidate = ID:\s*(\d+)/){
				my $cand_id=$1;
				my $conf;
				if(/Conf:\s*(\d(\.\d+)?)/){
					$conf=$1;
				};
				if(/Taxa:\s*\{\s*(.+?)\s*\}/){
					my @tmp=split(/\s+/, $1);
					my @taxa;
					if(defined $ra_rmtaxa){
						foreach my $str(@tmp){
							if(defined $removed_taxa{$str}){
								@taxa=();
								last;
							}else{
								push @taxa,$str;
							}
						}
					}else{@taxa=@tmp};
					next if @taxa==0;
					my @clades=find_clades($tree,\@taxa);
					foreach my $node(@clades){
						if(!defined $reassortments{$node->get_name}){
							$reassortments{$node->get_name}=new_reassortment($node);
						}else{
							my $n=$reassortments{$node->get_name}->nobs;
							if($seg_name ne $reassortments{$node->get_name}->info->[-1]->dataset){
								#new segment
								$reassortments{$node->get_name}->nobs($n+1);
							}
						};
						my $r_info=ReassortmentInfo->new();
						$r_info->id($cand_id);
						$r_info->dataset($seg_name);
						$r_info->conf($conf);
						$r_info->cardinality(scalar(@taxa));
						push @{$reassortments{$node->get_name}->info},$r_info;
						$reassortments{$node->get_name}->conf($conf) if($conf>$reassortments{$node->get_name}->conf);
					};
				};
			};
		};
		close REPORT;
	};
	foreach my $name(keys %reassortments){
		delete $reassortments{$name} if ($reassortments{$name}->nobs<$min_robs)&&($reassortments{$name}->conf<$min_rconf);
	};
	return %reassortments;
};

sub read_reassortment_catalog{
	my ($catalog_fn,$tree,$dataset_name,$ra_rmtaxa)=@_;
	my %removed_taxa;
	if(defined($ra_rmtaxa)&&(@{$ra_rmtaxa}>0)){
		foreach my $str(@{$ra_rmtaxa}){
			$removed_taxa{$str}=1;
		}
	}
	my %reassortments;
	open REPORT, "<$catalog_fn" or die "\nUnable to open input file: $catalog_fn!";
	my $cand_id=0;
	while(<REPORT>){
		my @line=split /:/;
		next if @line!=3;
		$cand_id++;
		for(my $i=0;$i<3;$i++){
			$line[$i]=~s/^\s+//;
			$line[$i]=~s/\s+$//;
		};
		my @tmp=split(/\s+/,$line[0]);
		my @taxa;
		if(defined $ra_rmtaxa){
			foreach my $str(@tmp){
				if(defined $removed_taxa{$str}){
					@taxa=();
					last;
				}else{
					push @taxa,$str;
				}
			}
		}else{@taxa=@tmp}
		next if @taxa==0;
		my @datasets=split(/\|/,$line[2]);
		for(my $i=0;$i<2;$i++){
			$datasets[$i]=~s/^\s+//;
			$datasets[$i]=~s/\s+$//;
		};
		if($datasets[0]=~/$dataset_name/){
			@datasets=split(/\s+/,$datasets[1]);
		}elsif($datasets[1]=~/$dataset_name/){
			@datasets=split(/\s+/,$datasets[0]);
		}else{next;};
		my @clades=find_clades($tree,\@taxa);
		foreach my $node(@clades){
			if(!defined $reassortments{$node->get_name}){
				$reassortments{$node->get_name}=new_reassortment($node);
			};
			foreach my $seg_name(@datasets){
				my $r_info=ReassortmentInfo->new();
				$r_info->id($cand_id);
				$r_info->dataset($seg_name);
				$r_info->cardinality(scalar(@taxa));
				push @{$reassortments{$node->get_name}->info},$r_info;
			};
		};
	};
	close REPORT;
	return %reassortments;
};

sub print_reassortment_info{
	my $fh=shift;
	my $rh_reassortments=shift;
	foreach my $name(keys %{$rh_reassortments}){
		print $fh "\n\tbegin clade;\n\t\tNode=\"$name\"\n\t\tGiraf={";
		my $r_info=$rh_reassortments->{$name}->info->[0];
		my $id=$r_info->id;
		my $dset=$r_info->dataset;
		print $fh "$dset:$id";
		for(my $i=1;$i<@{$rh_reassortments->{$name}->info};$i++){
			$r_info=$rh_reassortments->{$name}->info->[$i];
			$id=$r_info->id;
			$dset=$r_info->dataset;
			print ", $dset:$id";
		};
		my $conf=$rh_reassortments->{$name}->conf;
		print $fh "}\n\t\tConf=$conf\n\tend;"
	}
};

sub get_outgroup{
	my $tree=shift;
	my $outgrNode=$tree->get_root_node()->get_child(0);
	$outgrNode=$tree->get_root_node()->get_child(1) unless $outgrNode->is_terminal;
	$outgrNode=undef unless $outgrNode->is_terminal;
	return $outgrNode;
};

sub set_outgroup{
	my ($tree,$outgroup_str)=@_;
	my $outgrNode=get_outgroup($tree);
	my $node=$tree->find_node($outgroup_str);
	my $oldRoot=$tree->get_root_node();
#my $n1=@{$oldRoot->get_children};
	my $t=$node->set_root_below() unless(defined($outgrNode) && ($outgrNode == $node));
#my $n2=@{$oldRoot->get_children};
#print "\nn1=$n1\tn2=$n2";
	$oldRoot->collapse;
	return $t;
};


#Sets branch length in the '$tree' equal to branches in the another tree with relevant topology (identical or not contradictional) 
#from the NEWICK file '$ref_tree_fn'
sub reweight_tree{
	my ($tree,$ref_tree_fn, $outgroup,$ra_rmtaxa)=@_;
	my $treeio = Bio::TreeIO->new(-format => 'newick', -file => $ref_tree_fn);
	my $ref_tree   = $treeio->next_tree;
	if(!defined $outgroup){
		my $outgrNode=get_outgroup($tree);
		$outgroup=$outgrNode->get_name;
	};
	my $node=$ref_tree->find_node($outgroup);
	if(defined($node->ancestor->ancestor)){
		if($node->branch_length()){
			$ref_tree->reroot_at_midpoint($node);
		}else{
			$ref_tree->reroot($node->ancestor);
		}
	}
	$ref_tree=Bio::Phylo::Forest::Tree->new_from_bioperl($ref_tree);
	$ref_tree=$ref_tree->prune_tips($ra_rmtaxa) if defined $ra_rmtaxa && @{$ra_rmtaxa}>0;
	my %node_map;
	struct NodeInfo =>{
		xref => '$',
		#node =>'$',
		ntips => '$'
	};
	my %ref_clade_size;
	$ref_tree->visit_depth_first(
		'-pre_daughter' => sub{
			my $node=shift;
			$ref_clade_size{$node}=0;
		},
		'-in' => sub{ 
			my $node=shift;
			my $anc=$node->get_parent;
			if($anc){
				if($node->is_terminal){
					$ref_clade_size{$node}=1;
					my $node_inf=NodeInfo->new();
					$node_inf->ntips(1);
					$node_inf->xref($node);
					$node_map{$node->get_name}=$node_inf;
				};
				$ref_clade_size{$anc}+=$ref_clade_size{$node};
			}
		}
	);
	$tree->visit_depth_first(
	'-pre_daughter' => sub{
		my $node_inf=NodeInfo->new();
		$node_inf->ntips(0);
		my $node=shift;
		$node_map{$node->get_name}=$node_inf;
	},
	'-post_daughter' => sub{
		my $node=shift;
		my @xnode;
		foreach my $chnode(@{$node->get_children}){
			push @xnode, $node_map{$chnode->get_name}->xref;
		};
		my $xanc=$ref_tree->get_mrca(\@xnode);
		$node_map{$node->get_name}->xref($xanc);
	},
	'-in' => sub{ 
		my $node=shift;
		my $str=$node->get_name;
		my $anc=$node->get_parent;
		if($anc){
			my $ntips=$node_map{$anc->get_name}->ntips;
			$ntips+=$node_map{$str}->ntips;
			$node_map{$anc->get_name}->ntips($ntips);
		};
		my $xnode=$node_map{$str}->xref;
		if($node_map{$str}->ntips==$ref_clade_size{$xnode}){
			$node->set_branch_length($xnode->get_branch_length);
		}else{
			$node->set_branch_length(0);
		};
	}
	);
#Distance test
#my %taxa;
#my @tax_names;
#foreach my $node(@{$tree->get_terminals}){
#	$taxa{$node->get_name}=$node;
#	push @tax_names,$node->get_name;
#};
#my %taxa1;
#foreach my $node(@{$ref_tree->get_terminals}){
#	$taxa1{$node->get_name}=$node;
#};
#print "\nReweight tree distance test";
#	for(my $j=0;$j<@tax_names;$j++){
#		my $d=$taxa{$outgroup}->calc_patristic_distance($taxa{$tax_names[$j]});
#		my $d1=$taxa1{$outgroup}->calc_patristic_distance($taxa1{$tax_names[$j]});
#		print "\n$outgroup\t$tax_names[$j]\t$d\t$d1" if abs($d-$d1)>1e-5;
#	}
#print "\nDone!";
#exit;
#end distance test. Passed!
};

#This function splits a clade that contains taxons from predefined subset @{$ra_taxset}
#into clades which do not contain taxons which were not included in that subset.
sub split_clade{
	my ($mca,$ra_taxset)=@_;
	my %tax_names;
	my @splits;
	foreach my $tax(@{$ra_taxset}){
		$tax_names{$tax->get_name}=1;
	};
	my %leaf_counts;
	$mca->visit_depth_first(
		-pre_daughter => sub{
			my $node=shift;
			$leaf_counts{$node->get_name}=[(0,0)];
		},
		-no_daughter => sub{
			my $node=shift;
			$leaf_counts{$node->get_name}=[(1,0)];
			if(defined($tax_names{$node->get_name})){
				$leaf_counts{$node->get_name}->[1]=1;
			};
		},
		-in => sub{
			my $node=shift;
			my $pnode=$node->get_parent;
			if($pnode){
				$leaf_counts{$pnode->get_name}->[0]+=$leaf_counts{$node->get_name}->[0];
				$leaf_counts{$pnode->get_name}->[1]+=$leaf_counts{$node->get_name}->[1];
			}
		}
	);
	if($leaf_counts{$mca->get_name}->[0]!=$leaf_counts{$mca->get_name}->[1]){
		my @nodes=($mca);
		while(@nodes){
			my $node=shift @nodes;
			foreach my $chnode(@{$node->get_children}){
				if($leaf_counts{$chnode->get_name}->[0]==$leaf_counts{$chnode->get_name}->[1]){
					push @splits, $chnode;
				}else{
					push @nodes, $chnode;
				}
			}
		}
	};
	return @splits;
};	

#The function searches for each node in the tree towards the root direction the nearest ancestor from the node set '$rh_clade_set'
#and calculates distance.
#The optional parameter '$orh_bref' is a hash reference for returning the nearest ancestor from the '$rh_clade_set' for each branch in the tree '$tree'.
#The optional parameter '$dist_mode' defines the 5 posiible distance modes between branches in the rooted tree:
#	='BB' - begin - begin,
#	='BE' - begin - end,
#	='EB' - end - begin. The default mode,
#	='EE' - end - end,
#	='MB'- middle - begin,
#	='ME'- middle - end,
#	='BM'- begin - middle,
#	='MM'- middle - middle,
sub get_tree_markup_distribution{
	my ($tree,$rh_clade_set,$orh_bref,$dist_mode)=@_;
	my $inf_distance=100000;
	my %distances;
	my %bref;
	my $root=$tree->get_root_node;
	my @dist_mode=('E','B');
	if(defined $dist_mode){
		die "\nError in get_tree_markup_distribution(): unknown distance mode (the template is [BEM]{2})!" unless $dist_mode=~/[BEM]{2}/;
		@dist_mode=split "", $dist_mode;
	};
	$distances{$root->get_name}=$inf_distance;
	$tree->visit_breadth_first(
		-pre_daughter   => sub{
			my $node=shift;
			my $w=$distances{$node->get_name};
			foreach my $chnode(@{$node->get_children}){
				if(!defined($rh_clade_set->{$chnode->get_name})){
					if($w!=$inf_distance){
						$distances{$chnode->get_name}=$w+$chnode->get_branch_length;
						$bref{$chnode->get_name}=$bref{$node->get_name};
					}else{
						$distances{$chnode->get_name}=$inf_distance;
					};
				}else{
					if($dist_mode[0] eq 'B'){
						#from the edge beginning
						$distances{$chnode->get_name}=$chnode->get_branch_length;
					}elsif($dist_mode[0] eq 'M'){
						#from the edge middle
						$distances{$chnode->get_name}=$chnode->get_branch_length/2;
					}else{
						$distances{$chnode->get_name}=0;
					};
					$bref{$chnode->get_name}=$chnode;
				}
			}
		},
		-post => sub{
			my $node=shift;
			if($distances{$node->get_name}==$inf_distance){
				$distances{$node->get_name}='INF';
			}elsif(!defined($rh_clade_set->{$node->get_name})){
				if($dist_mode[1] eq 'B'){
					$distances{$node->get_name}-=$node->get_branch_length;
				}elsif($dist_mode[1] eq 'M'){
					$distances{$node->get_name}-=$node->get_branch_length/2;
				}
			}else{
				$distances{$node->get_name}=0;
			};
		}
	);
	%{$orh_bref}=%bref if(defined $orh_bref);
	return %distances;
};

#The function calculates distances for each node from the node set '$rh_clade_set' to the nearest ancestor from the same node set.
#The optional parameter '$orh_bref' is a hash reference for returning the nearest ancestor from the '$rh_clade_set' for each branch in the tree '$tree'.
#The optional parameter '$dist_mode' defines the 5 posiible distance modes between branches in the rooted tree:
#	='BB' - begin - begin,
#	='BE' - begin - end,
#	='EB' - end - begin. The default mode,
#	='EE' - end - end,
#	='MB'- middle - begin,
#	='ME'- middle - end,
#	='BM'- begin - middle,
#	='MM'- middle - middle,
sub get_tree_mrca_distances{
	my ($tree,$rh_clade_set,$orh_bref,$dist_mode)=@_;
	my $inf_distance=100000;
	my %distances;
	my %bref;
	my %rdist;
	my %rbref;
	my $root=$tree->get_root_node;
	my @dist_mode=('E','B');
	if(defined $dist_mode){
		die "\nError in get_tree_markup_distribution(): unknown distance mode (the template is [BEM]{2})!" unless $dist_mode=~/[BEM]{2}/;
		@dist_mode=split "", $dist_mode;
	};
	$distances{$root->get_name}=$inf_distance;
	$tree->visit_breadth_first(
		-pre_daughter   => sub{
			my $node=shift;
			my $w=$distances{$node->get_name};
			foreach my $chnode(@{$node->get_children}){
				if(!defined($rh_clade_set->{$chnode->get_name})){
					if($w!=$inf_distance){
						$distances{$chnode->get_name}=$w+$chnode->get_branch_length;
						$bref{$chnode->get_name}=$bref{$node->get_name};
					}else{
						$distances{$chnode->get_name}=$inf_distance;
					};
				}else{
					if($w!=$inf_distance){
						$rdist{$chnode->get_name}=$w+$chnode->get_branch_length;
						$rbref{$chnode->get_name}=$bref{$node->get_name};
					}else{
						$rdist{$chnode->get_name}=$inf_distance;
					};
					if($dist_mode[0] eq 'B'){
						#from the edge beginning
						$distances{$chnode->get_name}=$chnode->get_branch_length;
					}elsif($dist_mode[0] eq 'M'){
						#from the edge middle
						$distances{$chnode->get_name}=$chnode->get_branch_length/2;
					}else{
						$distances{$chnode->get_name}=0;
					};
					$bref{$chnode->get_name}=$chnode;
				}
			}
		},
		-post => sub{
			my $node=shift;
			if($distances{$node->get_name}==$inf_distance){
				$distances{$node->get_name}='INF';
			}elsif(!defined($rh_clade_set->{$node->get_name})){
				if($dist_mode[1] eq 'B'){
					$distances{$node->get_name}-=$node->get_branch_length;
				}elsif($dist_mode[1] eq 'M'){
					$distances{$node->get_name}-=$node->get_branch_length/2;
				}
			}else{
				$distances{$node->get_name}=0;
			}
			if(defined $rdist{$node->get_name}){
				if($rdist{$node->get_name}==$inf_distance){
					$rdist{$node->get_name}='INF';
				}else{
					if($dist_mode[1] eq 'B'){
						$rdist{$node->get_name}-=$node->get_branch_length;
					}elsif($dist_mode[1] eq 'M'){
						$rdist{$node->get_name}-=$node->get_branch_length/2;
					}
				}
			}
		}
	);
	%{$orh_bref}=%rbref if(defined $orh_bref);
	return %rdist;
};

sub get_dead_nodes{
	my ($tree,$ra_remove_taxids)=@_;
	my %remove_ids;
	foreach my $nodeid(@{$ra_remove_taxids}){
		$remove_ids{$nodeid}=1;
	};
	my @remove_nodes;
	foreach my $node(@{$tree->get_terminals}){
		push @remove_nodes,$node if(defined $remove_ids{$node->get_name});
	};
	my %parents;
	do{
		foreach my $node(@remove_nodes){
			my $pnode=$node->get_parent;
			$parents{$pnode->get_name}=[(0,$pnode)] unless defined $parents{$pnode->get_name};
			$parents{$pnode->get_name}->[0]--;
		};
		@remove_nodes=();
		foreach my $nodeid(keys %parents){
			my $node=$parents{$nodeid}->[1];
			my $n=$parents{$nodeid}->[0]+@{$node->get_children};
			if($n==0){
				push @remove_nodes, $node;
				delete $parents{$nodeid};
			};
			$remove_ids{$node->get_name}=1 if($n<2);
		}
	}while(@remove_nodes);
	return %remove_ids;
};

1;