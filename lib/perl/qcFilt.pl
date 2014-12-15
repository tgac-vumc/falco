#!/usr/bin/perl

use strict;

my $qcQual = shift;
my $output = shift;

my $minDp = shift || 100;
my $minVar = shift || 10;
my $minIns = shift || $minVar;
my $minDel = shift || $minVar;
my $minVarQ = shift || 100;
my $minInsQ = shift || 200;
my $minDelQ = shift || 200;

open QC, "<$qcQual";
my $head = readline QC;
chomp;
$head =~ s/\s+$//;
my @rhead = split(/\t/, $head);
my $coln = 0;
my %col = map { $_ => $coln++ } @rhead;
open OUT, ">$output.qc.ann.qual.filt.txt"; 
print OUT "$head\n";

while (<QC>) {
	my $line = $_;
	my @row = split(/\t/, $_);
	my $key = $row[$col{"X.chr"}] . ":" . $row[$col{pos}];
	next if ($row[$col{dp}] < $minDp);
	# Var
	if ($row[$col{qScore}] >= $minVarQ && $row[$col{nVar}] > $minVar) {
		print OUT $line;		
	}
	# Ins
	if ($row[$col{qScoreI}] >= $minInsQ) {
		$row[$col{ins}] =~ /:(\d+)/;
		if ($1 >= $minIns) {
			print OUT $line;
		}
	}
	# Del
	if ($row[$col{qScoreD}] >= $minDelQ) {
		$row[$col{del}] =~ /:(\d+)/;
		if ($1 >= $minDel) {
			print OUT $line;
		}
	}
}

close OUT;
