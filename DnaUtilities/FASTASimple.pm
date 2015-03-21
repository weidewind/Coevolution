package DnaUtilities::FASTASimple;
#Note! Do not use this module for huge files!

use strict;
use vars qw(@ISA @EXPORT_OK $VERSION);
use Exporter;
$VERSION = 1.00; # Or higher
@ISA = qw(Exporter);
#@EXPORT = qw(); # Symbols to autoexport (:DEFAULT tag)
@EXPORT_OK = qw(fetchFASTA_by_ids);

#returns a number of fetched sequences
sub fetchFASTA_by_ids{
	my ($seqfile,$rin_ids,$rout_snames,$rout_seqs)=@_;
	die "\n2-nd [in] parameter with ARRAYREF to sequense is list not defined!" unless $rin_ids;
	die "\n3-rd [out] parameter with ARRAYREF to sequence name list not defined!" unless $rout_snames;
	die "\n3-rd [out] parameter with ARRAYREF to sequence list not defined!" unless $rout_seqs;
	
	@{$rout_snames}=();
	@{$rout_seqs}=();
	my %ids  = ();
	foreach my $id(@{$rin_ids}){
		$ids{$id}+=1;
	};

	local $/ = "\n>";  # read by FASTA record

	open FASTA, "< $seqfile" or die "\nUnable to open input file: $seqfile!";
	my $n=0;
	while (<FASTA>) {
    chomp;
    my $seq = $_;
    my ($id) = $seq =~ /^>*([^\s]+)/;  # parse ID as first word in FASTA header
    if (exists($ids{$id})) {
    		$n++;
        $seq =~ s/^>*(.+)\n//;  # remove FASTA header
        push @{$rout_snames}, $1;
        $seq =~ s/\n//g;  # remove endlines
        $seq =~ s/\s//g;	# remove white spaces
        push @{$rout_seqs}, $seq;
    }
	}
	close FASTA;
	$n;
};
