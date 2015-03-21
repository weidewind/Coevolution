package DnaUtilities::iupac;

use strict;
use vars qw(@ISA @EXPORT $VERSION);
use Exporter;
$VERSION = 1.00; # Or higher
@ISA = qw(Exporter);
@EXPORT = qw(ch2code ch2idx idx2ch code2iupac base2iupac iupac_join iupac_union iupac_vol base_cmp); # Symbols to autoexport (:DEFAULT tag)

#Utilities
#Obsolete
sub ch2code{
	my $s=$_[0];
	if($s =~/[Aa]/){
		return 0;
	}elsif($s =~/[Cc]/){
		return 1;
	}elsif($s =~/[Gg]/){
	  return 2;
	}elsif($s =~/[TtUu]/){
		return 3;
	};	
	return -1;
};

sub ch2idx{
	return ch2code @_;
};

sub idx2ch{
	my $idx=$_[0];
	if($idx=~m/^\s*([0123])\s*$/i){
		$idx=$1;
	}else{
		die "\nWrong base index: $idx";
	};
	my @a=('a','c','g','t');
	return $a[$idx];
};
 
#Private
sub code2iupac{
	my $code=$_[0];
	if($code eq "1000"){
		return 'a';
	}elsif($code eq "0100"){
		return 'c';
	}elsif($code eq "0010"){
		return 'g';
	}elsif($code eq "0001"){
		return 't';
	}elsif($code eq "1100"){
		return 'm';
	}elsif($code eq "1010"){
		return 'r';
	}elsif($code eq "1001"){
		return 'w';
	}elsif($code eq "0110"){
		return 's';
	}elsif($code eq "0101"){
		return 'y';
	}elsif($code eq "0011"){
		return 'k';
	}elsif($code eq "1110"){
		return 'v';
	}elsif($code eq "0111"){
		return 'b';
	}elsif($code eq "1011"){
		return 'd';
	}elsif($code eq "1101"){
		return 'h';
	}elsif($code eq "1111"){
		return 'n';
	};
	return '-';
};

sub base2iupac{
	my $sym=shift;
	$sym=lc $sym;
	my @r=();
	if($sym eq 'a'){
		@r=(1,0,0,0);
	}elsif($sym eq 'c'){
		@r=(0,1,0,0);
	}elsif($sym eq 'g'){
		@r=(0,0,1,0);
	}elsif($sym eq 't'){
		@r=(0,0,0,1); ;
	}elsif($sym eq 'm'){
		@r=(1,1,0,0);
	}elsif($sym eq 'r'){
		@r= (1,0,1,0);
	}elsif($sym eq 'w'){
		@r= (1,0,0,1);
	}elsif($sym eq 's'){
		@r= (0,1,1,0);
	}elsif($sym eq 'y'){
		@r= (0,1,0,1);
	}elsif($sym eq 'k'){
		@r= (0,0,1,1);
	}elsif($sym eq 'v'){
		@r= (1,1,1,0);
	}elsif($sym eq 'b'){
		@r= (0,1,1,1);
	}elsif($sym eq 'd'){
		@r= (1,0,1,1);
	}elsif($sym eq 'h'){
		@r= (1,1,0,1);
	}elsif($sym eq 'n'){
		@r= (1,1,1,1);
	}elsif($sym eq '-'){
		@r=(0,0,0,0);
	};
	if (wantarray()) {
		return @r;
	}elsif (defined wantarray()) {
		my $code=join "",@r;
		return $code;
	}else {
		# void context
		print "\nWarning BCVIndel::Patterns::base2iupac called in void context!";
		return;
	};
};

sub iupac_join{
	my ($r1,$r2)=@_;
	return ($r1->[0]&&$r2->[0],$r1->[1]&&$r2->[1],$r1->[2]&&$r2->[2],$r1->[3]&&$r2->[3]);
};

sub iupac_union{
	my ($r1,$r2)=@_;
	return ($r1->[0]||$r2->[0],$r1->[1]||$r2->[1],$r1->[2]||$r2->[2],$r1->[3]||$r2->[3]);
};

sub iupac_vol{
	my $r=shift;
	return $r->[0]+$r->[1]+$r->[2]+$r->[3];
};

sub base_cmp{
	my ($a,$b)=@_;
	my @A=base2iupac($a);
	my @B=base2iupac($b);
	my @I=iupac_join(\@A,\@B);
	my @J=iupac_union(\@A,\@B);
	return iupac_vol(\@I)/iupac_vol(\@J);
};