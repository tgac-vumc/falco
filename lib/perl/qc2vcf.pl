#!/usr/bin/perl

use strict;
use File::Basename;

my $qc = shift;
my $base = shift || undef;

if ($base) {
	$base = basename($qc, "qc.ann.qual.filt.txt")
}

my $minDp = 100;
my $minVarDp = 5;
my $minVaf = 0.01;

my $header = qq(##fileformat=VCFv4.1
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Depth">
##FORMAT=<ID=AD,Number=2,Type=Integer,Description="Allele Distribution">
##FORMAT=<ID=AB,Number=1,Type=Float,Description="Allele Ballance nref/nref+nvar">
##INFO=<ID=TARGET,Number=1,Type=String,Description="TSACP Target name">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	$base
);

my $chr = 0; 
my $pos = 1;
my $target = 2; 
my $dp = 3;
my $ref = 4;
my $var = 5;
my $nref = 6;
my $nvar = 7;
my $nA = 8;
my $nC = 9;
my $nG = 10;
my $nT = 11;
my $nN = 12;
my $ins = 13;
my $del = 14;
my $pval = 15;
my $pvala = 16;
my $varQual = 17;
my $insQual = 18;
my $delQual = 19;

my %acgt = (A => 0, C => 1, G => 2, T => 3);
print $header;
open QC, "<$qc";
while (<QC>) {
	chomp;
	next if (/^#/);
	next if (/^X\./);
	my @row = split(/\t/, $_);
	next if ($row[$target] eq ".");

	my $VAR = substr($row[$var], 0, 1);
	my $REF = $row[$ref];
	my $nvarall = $row[$nvar];
	my $indel = $row[$ins] . $row[$del];	
	if ($row[$nvar] > 0 || $indel eq "..") {
		my %ord = ( A => $row[$nA], C => $row[$nC], G => $row[$nG], T => $row[$nT], N => $row[$nN] );
		$ord{uc($REF)} = -1; # Exclude ref from sort
		$ord{N} = -1; # Exclude N!
		my @order = sort { $ord{$b} <=> $ord{$a}  } keys(%ord);
		my $dpVar = $ord{$order[0]};
		my $VAR = $order[0];
		my $skip = 0;
		my $tot = $row[$nref] + $dpVar;
		my $vaf = ($tot > 0)?$dpVar / $tot:"NA"; # Excluding other alleles...
		#$skip++ if ($tot == 0);
		#$skip++ if ($ord{$order[0]} == $ord{$order[1]}); # Multiallelic, skipping for now
		#$skip++ if ($dpVar == 0); # Skip if dpVar is zero, can occur when N is the only variant. Above line also protects for this.
		#$skip++ if ($dpVar < $minVarDp);
		#$skip++ if ($row[$dp] < $minDp);
		#$skip++ if ($vaf < $minVaf);
		my $gt = "$row[$ref]/$VAR";
		my $qual = $row[$varQual] || 0;
		$skip++ if ($qual == 0);
		my $filter = 0;
		my $info = "TARGET=$row[$target]";
		my $format = "DP:AD:AB";
		my $sample = "$row[$dp]:$row[$nref],$dpVar:$vaf";
		print join("\t", @row[$chr, $pos], ".", $REF, $VAR, $qual, $filter, $info, $format, $sample) . "\n" if ($skip == 0);
	}
	# INDEL
	if ($row[$ins] ne ".") {
		my $skip = 0;
		# Modify coordinte to conform to VCF?
		my %inss = split(/[\|\:]/, $row[$ins]);
		my @keys = sort { $inss{$b} <=> $inss{$a} } keys %inss;
#		if ($inss{$keys[0]} > $row[$nvar]) {
#			$nvarall = $inss{$keys[0]};
#			$VAR = $keys[0];
#			$row[$dp] += $nvarall;
#		}
		my $nins = $inss{$keys[0]};
		my $DP = $row[$dp];
		my $nREF = $DP - $nins; 
		#$skip++ if ($nins < $minVarDp);
		#$skip++ if ($DP < $minDp);
		my $vaf = $nins / $row[$dp];
		#$skip++ if ($vaf < $minVaf);
		my $gt = "$row[$ref]/$keys[0]";
		my $qual = $row[$insQual];
		$skip++ if ($qual == 0);
		my $filter = 0;
		my $info = "TARGET=$row[$target]";
		my $format = "DP:AD:AB";
		my $sample = "$DP:$nREF,$nins:$vaf";
		print join("\t", @row[$chr, $pos], ".", $REF, $keys[0], $qual, $filter, $info, $format, $sample) . "\n" if ($skip == 0);
	}
	if ($row[$del] ne ".") {
		my $skip = 0;
		# Modify coordinte to conform to VCF?
		my %dels = split(/[\|\:]/, $row[$del]);
		my @keys = sort { $dels{$b} <=> $dels{$a} } keys %dels;
#		if ($dels{$keys[0]} > $row[$nvar]) {
#			$nvarall = $dels{$keys[0]};
#			$REF = $keys[0];
#			$row[$dp] += $nvarall;
#		}
		my $ndel = $dels{$keys[0]};
		my $DP = $row[$dp];
		my $nREF = $DP - $ndel;
		#$skip++ if ($ndel < $minVarDp);
		#$skip++ if ($DP < $minDp);
		my $vaf = $ndel / $DP;
		#$skip++ if ($vaf < $minVaf);
		my $gt = "$keys[0]/$row[$ref]";
		my $qual = $row[$delQual];
		$skip++ if ($qual == 0);
		my $filter = 0;
		my $info = "TARGET=$row[$target]";
		my $format = "DP:AD:AB";
		my $sample = "$DP:$nREF,$ndel:$vaf";
		print join("\t", @row[$chr, $pos], ".", $keys[0], $row[$ref], $qual, $filter, $info, $format, $sample) . "\n" if ($skip == 0);
	}
}
close QC;

