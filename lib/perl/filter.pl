#!/usr/bin/perl -w

use strict;

my $file = shift;
my $canonicals = shift;
my $clinvar = shift;
my $cosmic = shift;
my $cosmicNC = shift;

my $minDP = 100;
my $minVAR = 10;
my $minVAF = 0.01;

my %can = ();
my %nm = ();
my %cvHsh = ();

print STDERR "Reading in canonicals,...";
open CAN, "<$canonicals";
while (<CAN>) {
	chomp;
	my @row = split(/\t/, $_);
	$can{$row[0]} = @row[1 .. $#row];
	my $pri = 0;
	$nm{$_} = $pri++ foreach @row[1 .. $#row];
}
close CAN;

print STDERR "Reading clinvars...\n";
open CV, "<$clinvar";
while (<CV>) {  last if (/#CHROM/); }
while (<CV>) {
        chomp;
        my @row = split(/\t/, $_);
        push @{$cvHsh{"chr".$row[0].":".$row[1]}}, [@row];
}
close CV;

print STDERR "Reading COSMIC...\n";
open CV, "<$cosmic";
while (<CV>) {  last if (/#CHROM/); }
while (<CV>) {
        chomp;
        my @row = split(/\t/, $_);
        push @{$cvHsh{"chr".$row[0].":".$row[1]}}, [@row];
}
close CV;

print STDERR "Reading COSMICNC...\n";
open CV, "<$cosmicNC";
while (<CV>) {  last if (/#CHROM/); }
while (<CV>) {
        chomp;
        my @row = split(/\t/, $_);
        push @{$cvHsh{"chr".$row[0].":".$row[1]}}, [@row];
}

close CV;
print STDERR "done!\n";

print STDERR "Processing $file\n";

open F, "<$file";
my $head = readline(F);
chomp $head;

print $head . "\tvaf\n";

my @rowHead = split(/\t/, $head);
my $colN = 0;
my %col = map {
	$_ => $colN++ } @rowHead;

while (<F>) {
	chomp;
	my $line = $_;
	my @row = split(/\t/, $_);
	next unless (exists($can{$row[$col{Gene_Name}]}));
	my $trans = $row[$col{Transcript_ID}];
	$trans =~ s/\..*$//;
	
	if (not exists($nm{$trans})) {
		print STDERR "Unknown transcript: $row[$col{Gene_Name}] $trans\n";
		next;
	}
	elsif ($nm{$trans} != 0) {
		print STDERR "Non cannonical: $row[$col{Gene_Name}] $trans\n";
		next;
	}

	# Annotate clinincal variants and cosmic mutations

	my $chr = $row[$col{"#CHROM"}];
	my $pos = $row[$col{POS}];

	if (exists($cvHsh{"$chr:$pos"})) {
		my $rows = $cvHsh{"$chr:$pos"};
		my $cv = join("|", map { $_->[2] } @$rows); 
		$row[2] .= "|$cv";
	}		
		
	# Filter INTERGENIC
	next if ($row[$col{Context}] eq "INTERGENIC");
	
	# Filter SYNONYMOUS_CODING 
	# Dubbel mutatie BRAF 
#	next if ($row[$col{Context}] eq "SYNONYMOUS_CODING");
	
	# Filter INTRON
	next if ($row[$col{Context}] eq "INTRON");

	# Filter dbSNP 
#	next if ($row[2] ne ".");

	# Filter depth
	next if ($row[$col{DP}] < $minDP);
	
	my ($nref, $nvar) = split(/,/, $row[$col{AD}]);	
	# Filter var depth
#	print STDERR "DEBUG: " . $nvar . "\n"; sleep 1;
	next if ($nvar < $minVAR);
	
	my $vaf = $nvar / ($nref + $nvar);
#	print STDERR "DEBUG: " . $vaf . "\n";
	# Filter vaf
	next if ($vaf < $minVAF);

	print join("\t", @row, $vaf) . "\n";	
}
close F;
