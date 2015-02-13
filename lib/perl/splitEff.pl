#!/usr/bin/perl -w
#
use strict;
my $file = shift;
my $tag = shift || "NoTag";

open F, "<$file";
my $head = undef;
my %effIdx = ();
$effIdx{context} = 0;
while (<F>) {
	if (/.*ID=EFF.*(\(.*\)).*/) {
		my $str = $1;
#		print STDERR $str . "\n";
		$str =~ s/\[.*\]//g;
#		print STDERR $str . "\n";
		$str =~ s/^\(\s*(.*?)\s*\)$/$1/;
#		print STDERR $str . "\n";
		my $c = 1;
		$effIdx{lc($_)} = $c++ foreach split(/\s*\|\s*/, $str);
#		print STDERR join("\n", map { $_ . ":" . $effIdx{$_} } keys(%effIdx)) . "\n"; 
#		exit;
	}
	if (/#CHROM/) {
		$head = $_;
		last;
	}

}

chomp $head;
my @hrow = split(/\t/, $head);
my $colN = 0;
#my %col = map { print STDERR "$_ => $colN\n";$_ => $colN++ } @hrow;
my %col = map { $_ => $colN++ } @hrow;
# Add AD column
# Convert DP4 for samtools file
my $ADName = (grep { /\_AD/ } keys(%col))[0];
if (exists($col{DP4})) {
	$ADName = "DP4";
}
my $ADCol = $col{$ADName};

my $DPName = (grep { /\_DP/ } keys(%col))[0];
if (exists($col{DP})) {
	$DPName = "DP";
}
my $DPCol = $col{$DPName};

#print STDERR $ADName . "\t" . $col{$ADName} . " $ADCol\n";

#Effect_Impact | Functional_Class | Codon_Change | Amino_Acid_change| Amino_Acid_length | Gene_Name | Gene_BioType | Coding | Transcript | Exon [ | ERRORS | WARNINGS ] )
my @Col = qw/DP AD Context Effect_Impact Functional_Class Codon_Change Amino_Acid_change Amino_Acid_length Gene_Name Transcript_BioType Gene_Coding Transcript_ID Exon_Rank Tag/;

my @effCol = map {lc} qw/Context Effect_Impact Functional_Class Codon_Change Amino_Acid_change Amino_Acid_length Gene_Name transcript_BioType Gene_Coding Transcript_ID Exon_Rank/;
my @effColOrd = map {$effIdx{$_}} @effCol;

print join("\t", @hrow[0 ..5], @Col) . "\n"; #map { "eff$_" } 1 .. 10) . "\n";
while (<F>) {
	chomp;
	my @row = split(/\t/, $_);

	if ($ADName eq "DP4") {
		my @dp4row = split(/,/, $row[$ADCol]);
#		print STDERR join(":", @dp4row) . "\n"; sleep 1;
		$row[$ADCol] = ($dp4row[0] + $dp4row[1]) . "," . ($dp4row[2] + $dp4row[3]);		
	}	

#	print $row[$col{EFF}] . "\n";
	my @effs = split(/,/, $row[$col{EFF}]);
	foreach my $eff (@effs) {
		$eff =~ s/\)$//;
		my @effrow = split(/[\|\(\)]/, $eff, -1);
			
#		print STDERR scalar(@effrow) . " $eff\n";
#		print STDERR join(":", @effrow) . " $eff\n"; sleep 1;
		print join("\t", @row[0 .. 5, $DPCol, $ADCol], @effrow[@effColOrd], $tag) . "\n";

	}
}

close F;
