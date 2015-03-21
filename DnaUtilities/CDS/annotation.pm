package DnaUtilities::CDS::annotation;
#Note! Do not use this module for huge files!

use strict;
use vars qw(@ISA @EXPORT @EXPORT_OK $VERSION);
use Exporter;
use DnaUtilities::FASTASimple;
$VERSION = 1.00; # Or higher
@ISA = qw(Exporter);
@EXPORT = qw(alnpos2annotation read_annotation_xml); # Symbols to autoexport (:DEFAULT tag)
#@EXPORT_OK = qw();

use Class::Struct;
use XML::Simple;

#Annotation info
struct FeatureInfo => {
	name => '$',
	loci => '@'
};
struct AnnotationInfo => {
	skip_cols => '@',
	refseq_cols => '@',
	offset => '$',
	rfeatures => '@'
};

#This function returns the feature name and relative coordinates for an alignment position
sub alnpos2annotation{
	my ($segm_id,$apos,$ra_annotation)=@_;
	my @annotation=@{$ra_annotation};
	return $apos unless scalar @annotation;
	die "\nError alnpos2annotation(): wrong segment id: $segm_id!" unless defined $annotation[$segm_id];
	#Here is assumed that $apos is an aminoacid position in the spliced protein
	$apos*=3;
	my $pos=$apos;
	my $is_ins=0;
	if(@{$annotation[$segm_id]->skip_cols}){
		#restore initial alignment coordinates (this operation is reverse to the splicing operation)
		foreach my $indel(@{$annotation[$segm_id]->skip_cols}){
			if($apos>=$indel->[0]){
				$pos+=$indel->[1];
			}else{last;};
		};
		$apos=$pos;
	};
	#converts coordinates on alignment into the reference sequence coordinates
	if(@{$annotation[$segm_id]->refseq_cols}){
		$pos=0;
		foreach my $interval(@{$annotation[$segm_id]->refseq_cols}){
			if($apos>$interval->[1]){
				$pos+=$interval->[1]-$interval->[0]+1;
			}else{
				if($apos>=$interval->[0]){
					$pos+=$apos-$interval->[0]+1;
				}else{
					$is_ins=1;
				};
				last;
			};
		};
	};
	$pos+=$annotation[$segm_id]->offset;
	$pos-- if $annotation[$segm_id]->offset>0;
	my @features;
	#find features which contain the inspecting site
	if(@{$annotation[$segm_id]->rfeatures}){
		$apos=$pos;
		foreach my $rft(@{$annotation[$segm_id]->rfeatures}){
			$pos=0;
			next if $apos<$rft->loci(0)->[0]||$apos>$rft->loci(-1)->[1];
			foreach my $locus(@{$rft->loci}){
				if($apos>$locus->[1]){
					$pos+=$locus->[1]-$locus->[0]+1;
				}else{
					if($apos>=$locus->[0]){
						$pos+=$apos-$locus->[0]+1;
						push @features, [$rft->name(),$pos];
					};
					last;
				}
			}
		}
		if(@features){ 
			@features=sort {$a->[1]<=>$b->[1]} @features;
			$pos=$features[0]->[1];
			$pos/=3; #only ORF features are assumed
			return ($pos,$features[0]->[0],$is_ins);
		}else{
			$pos=$apos/3; #ORF is assumed
		};
	}else{
		$pos/=3; #ORF is assumed
	};
	return ($pos,"ORF",$is_ins);
}

#Read the annotation from the XML file
#Reading the XML annotation file
sub read_annotation_xml{
	my $annot_fn=shift;
	my @annotation;
	my $bcv_prj = XML::Simple->new();
	
	my $data   = $bcv_prj->XMLin($annot_fn, 
		ForceArray => ['exon','gene'],
		SuppressEmpty => 1,
		ContentKey => '-content' );
	#use Data::Dumper;
	#print Dumper($data);
	#exit;
	for(my $i=1;$i<=2;$i++){
		my $seg_info=$data->{Segment}->{$i};
		if(defined $seg_info){
			my $p=AnnotationInfo->new();
			if(defined $seg_info->{SkippedColumnsFile}){
				my $cols_fn=$seg_info->{SkippedColumnsFile};
				open INF, $cols_fn or die "\nUnable to open input file $cols_fn!";
				while(<INF>){
					push @{$p->skip_cols}, [$1,$2] if(/(\d+)\t(\d+)/);
				};
				close INF;
			};
			$p->offset($seg_info->{Offset});
			my $refseq_name=$seg_info->{RefSeqName};
			my $align_fn=$seg_info->{AlignmentFile};
			if(defined($refseq_name)&&defined($align_fn)){
				my @tmp;
				my @seqs;
				my @rfs=($refseq_name);
				if(!DnaUtilities::FASTASimple::fetchFASTA_by_ids($align_fn,\@rfs,\@tmp,\@seqs)){
					die "\nThe file $align_fn doesn't contain the sequence with id=$refseq_name!";
				};
				@tmp=split "", $seqs[0];
				my @I=(0,0);
				for(my $j=0;$j<@tmp;$j++){
					if($tmp[$j] ne "-"){
						$I[0]=$j+1 if $I[0]==0;
						$I[1]=$j+1;
					}elsif($I[1]){
						push @{$p->refseq_cols},[@I];
						@I=(0,0);
					};
				};
				if($I[1]){
					push @{$p->refseq_cols},[@I];
				};
			};
			if(defined $seg_info->{gene}){
				foreach my $gname(keys %{$seg_info->{gene}}){
					my $feature=FeatureInfo->new();
					$feature->name($gname);
					foreach my $exon(@{$seg_info->{gene}->{$gname}->{exon}}){
						my @I=($exon->{from},$exon->{to});
						push @{$feature->loci},[@I];
					};
					push @{$p->rfeatures}, $feature;
				};
			};
			$annotation[$i-1]=$p;
		};
	};	
	return @annotation;
}

1;