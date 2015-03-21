package DnaUtilities::GrouppedMSA;

use strict;
use Bio::SimpleAlign;
use Bio::LocatableSeq;
use Class::Struct;

#This class implements a Groupped Multiple Alignment

sub new{
	my $class=shift;
	my $self={};
	bless $self, $class;
	
	$self->{R_MSA}=undef;
	$self->{SEQ_NAMES}=[];
	$self->{SEQ_WEIGHTS}=[];
	$self->{GROUP_EDGES}=[];
	$self->{GROUP_NAMES}=[];
	$self->{GROUP_WEIGHTS}=[];
	$self->{LEADING_GAP_SYMBOL}="-";
	$self->{TRAILING_GAP_SYMBOL}="-";

	return $self;
};
	

sub readGFASTA{
	my $self=shift;
	my $fname=shift;
	my %myargs=(
		LEADING_GAP_SYMBOL => '-',
		TRAILING_GAP_SYMBOL => '-',
		NO_ALIGNMENT => 0, #read group and sequence name only
		@_
	);
	if($myargs{LEADING_GAP_SYMBOL} ne '-'){
		my $sym=split "", $myargs{LEADING_GAP_SYMBOL};
		if($sym ne "" && (!($sym =~/^s+$/))){
			$self->{LEADING_GAP_SYMBOL}=$sym;
		}else{
			$self->{LEADING_GAP_SYMBOL}=='-';
		};
	};
	if($myargs{TRAILING_GAP_SYMBOL} ne '-'){
		my $sym=split "", $myargs{TRAILING_GAP_SYMBOL};
		if($sym ne "" && (!($sym =~/^s+$/))){
			$self->{TRAILING_GAP_SYMBOL}=$sym;
		}else{
			$self->{TRAILING_GAP_SYMBOL}=='-';
		};
	};
	if(!$myargs{NO_ALIGNMENT}){
		$self->{R_MSA}=Bio::SimpleAlign->new();
	};
	open ALMN, "< $fname" or die "\nUnable to open input file: $fname";
	my $j=0;
	my $i=0;
	my $name=undef;
	my $str="";
	my $sstr="\\w*?.-";
	$sstr.=$self->{LEADING_GAP_SYMBOL} unless($sstr=~/$self->{LEADING_GAP_SYMBOL}/);
	$sstr.=$self->{TRAILING_GAP_SYMBOL} unless($sstr=~/$self->{TRAILING_GAP_SYMBOL}/);
	$sstr="[^".$sstr."]";
	while(<ALMN>){     
		if(/^\s*=/ || /^\s*>/){
			if($str ne ""){    
				#add sequence what have already been read into alignment 
				if($name=~/(0\.\d+)/){
					$self->{SEQ_WEIGHTS}->[$i]=$1;	
				}else{
					$self->{SEQ_WEIGHTS}->[$i]=0;
				};
				if($j>0){
					$self->{GROUP_EDGES}->[$j-1]->[0]=$i if($self->{GROUP_EDGES}->[$j-1]->[0]==-1);
					$self->{GROUP_EDGES}->[$j-1]->[1]=$i;
					$self->{GROUP_WEIGHTS}->[$j-1]+=$self->{SEQ_WEIGHTS}->[$i];
				};
				$self->{SEQ_NAMES}->[$i]=$name;
				if($self->{R_MSA}){
					if($self->{LEADING_GAP_SYMBOL} ne '-'){
						$str=~s/^($self->{LEADING_GAP_SYMBOL}+)//;
						my $i=length $1;
						if($i){
							my $lgap="";
							while($i){$lgap.="-";$i--};
							$str=$lgap.$str;
						};
					};
					if($self->{TRAILING_GAP_SYMBOL} ne '-'){
						$str=~s/($self->{TRAILING_GAP_SYMBOL}+)$//;
						my $i=length $1;
						if($i){
							my $lgap="";
							while($i){$lgap.="-";$i--};
							$str.=$lgap;
						};
					};
					my $l=length($str);
					while($str=~m/(-+)/g){
						$l-=length $1;
					};
					my $seq=Bio::LocatableSeq->new(-seq =>$str,
						-id => $i+1,
						-start => 1,
						-end => $l);
					$self->{R_MSA}->add_seq($seq);   
				};  
			};
			$str="";
			$name=undef;
			$i++;
		};
		if(/^\s*=\s*(.+?)\s*$/){
			#read name of new group
			$name=$1;
			push @{$self->{GROUP_EDGES}}, [(-1,-1)];
			push @{$self->{GROUP_NAMES}},$name;
			push @{$self->{GROUP_WEIGHTS}},0;
			$name="";
			$j++;
		}elsif(/^\s*>\s*(.+)\s*$/){
		  $name=$1;   
		  $name=~s/[^a-z\d]+\s*$//i;
		}else{
			s/$sstr//g;
			$str.=$_;
		};	
	}; 
	
	close ALMN;
	
	if($str ne ""){    
		#add sequence what have already been read into alignment 
		if($name=~/(0\.\d+)/){
			$self->{SEQ_WEIGHTS}->[$i]=$1;	
		}else{
			$self->{SEQ_WEIGHTS}->[$i]=0;
		};
		if($j>0){
			$self->{GROUP_EDGES}->[$j-1]->[0]=$i if($self->{GROUP_EDGES}->[$j-1]->[0]==-1);
			$self->{GROUP_EDGES}->[$j-1]->[1]=$i;
			$self->{GROUP_WEIGHTS}->[$j-1]+=$self->{SEQ_WEIGHTS}->[$i];
		};
		$self->{SEQ_NAMES}->[$i]=$name;     
		if($self->{R_MSA}){ 
			if($self->{LEADING_GAP_SYMBOL} ne '-'){
				$str=~s/^($self->{LEADING_GAP_SYMBOL}+)//;
				my $i=length $1;
				if($i){
					my $lgap="";
					while($i){$lgap.="-";$i--};
					$str=$lgap.$str;
				};
			};
			if($self->{TRAILING_GAP_SYMBOL} ne '-'){
				$str=~s/($self->{TRAILING_GAP_SYMBOL}+)$//;
				my $i=length $1;
				if($i){
					my $lgap="";
					while($i){$lgap.="-";$i--};
					$str.=$lgap;
				};
			};
			my $l=length($str);
			while($str=~m/(-+)/g){
				$l-=length $1;
			};
			my $seq=Bio::LocatableSeq->new(-seq =>$str,
				-id => $i+1,
				-start => 1,
				-end => $l);
			$self->{R_MSA}->add_seq($seq);  
		};
		$i++;
	};
	return ($j,$i);
};

sub _print{
	my $self=shift;
	for(my $j=0;$j<@{$self->{GROUP_NAMES}};$j++){
		print "=$self->{GROUP_NAMES}->[$j]\n";
		for(my $i=$self->{GROUP_EDGES}->[$j]->[0];$i<=$self->{GROUP_EDGES}->[$j]->[1];$i++){
			my $seq="";
			my $str="";
			if($self->{R_MSA}){
				$seq=$self->{R_MSA}->get_seq_by_pos($i+1);
				$str=$seq->seq();
			};
			my $sname=$self->{SEQ_NAMES}->[$i];
			print ">$sname\n$str\n";
		};
	};
};

sub no_alignment{
	my $self=shift;
	return 1 unless $self->{R_MSA};
	return 0;
};

sub get_group_number{
	my $self=shift;
	return scalar @{$self->{GROUP_NAMES}};
};

sub get_group_name{
	my $self=shift;
	my $grid=shift;
	$grid--;
	return undef unless(defined($grid) && $grid >=0 && $grid<@{$self->{GROUP_NAMES}});
	return $self->{GROUP_NAMES}->[$grid];
};

sub get_group_weight{
	my $self=shift;
	my $grid=shift;
	$grid--;
	return undef unless(defined($grid) && $grid >=0 && $grid<@{$self->{GROUP_NAMES}});
	return $self->{GROUP_WEIGHTS}->[$grid];
};

sub get_group_size{
	my $self=shift;
	my $grid=shift;
	$grid--;
	return undef unless(defined($grid) && $grid >=0 && $grid<@{$self->{GROUP_NAMES}});
	return 0 if $self->{GROUP_EDGES}->[$grid]->[1]==-1;
	return $self->{GROUP_EDGES}->[$grid]->[1]-$self->{GROUP_EDGES}->[$grid]->[0]+1;
};

struct SequenceInfo => {
	LocatableSeq => '$',
	SeqName => '$',
	SeqWeight => '$'
}; 

sub get_sequence_from_group{
	my $self=shift;
	my $grid=shift;
	my $seqid=shift;
	my $n=$self->get_group_size($grid);
	$grid--;
	return undef unless (defined($seqid) && $seqid>0 && $seqid<=$n);

	$seqid+=$self->{GROUP_EDGES}->[$grid]->[0];
	my $p=SequenceInfo->new();
	if($self->{R_MSA}){
		$p->LocatableSeq($self->{R_MSA}->get_seq_by_pos($seqid));
	};
	$p->SeqName($self->{SEQ_NAMES}->[$seqid-1]);
	$p->SeqWeight($self->{SEQ_WEIGHTS}->[$seqid-1]);
  return $p;
};

sub get_sequence_by_id{
	my $self=shift;
	my $seqid=shift;
	return undef unless (!$seqid && $seqid>0 && $seqid<=@{$self->{SEQ_NAMES}});
	my $p=SequenceInfo->new();
	if($self->{R_MSA}){
		$p->LocatableSeq($self->{R_MSA}->get_seq_by_pos($seqid));
	};
	$p->SeqName($self->{SEQ_NAMES}->[$seqid-1]);
	$p->SeqWeight($self->{SEQ_WEIGHTS}->[$seqid-1]);
  return $p; 
};

sub get_seqname_from_group{
	my $self=shift;
	my $grid=shift;
	my $seqid=shift;
	my $n=$self->get_group_size($grid);
	$grid--;
	return undef unless (defined($seqid) && $seqid>0 && $seqid<=$n);

	$seqid+=$self->{GROUP_EDGES}->[$grid]->[0];
	return $self->{SEQ_NAMES}->[$seqid-1];
};

sub get_seqname_by_id{
	my $self=shift;
	my $seqid=shift;
	return undef unless (!$seqid && $seqid>0 && $seqid<=@{$self->{SEQ_NAMES}});
	return $self->{SEQ_NAMES}->[$seqid-1];
};


1