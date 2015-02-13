#!/usr/bin/perl -w

use strict;

my @vcfs = @ARGV;

my @info = ();
my @format = ();
my $header = "";

foreach my $vcf (0 .. $#vcfs) {
	open V, "<$vcfs[$vcf]";
	while (<V>) {
		if (/^#CHROM/) {
			if ($vcf == 0) {
				chomp;
				my @row = split(/\t/, $_);
				print $header;
				print join("\t", @row[0 .. 6], @info) . "\t";
				my @sampleRow = ();
				foreach my $s (9 .. $#row) {
#					print join("\t", map { "S$s" . "_$_"  } @format) . "\n";
					push @sampleRow, map { $row[$s] . "_$_" } @format; 		
				}
				print join("\t", @sampleRow) . "\n";
			}
			last;
		}
		if (/^##INFO/ && ($vcf == 0)) {
			/ID=(.*?),/;
			push @info, $1;	
		}
		elsif (/^##FORMAT/ && ($vcf == 0)) {
			/ID=(.*?),/;
			push @format, $1;
		}
		$header .= $_;
	}

	while (<V>) {
		print join("\t", @{ &readTags(\$_) }) . "\n";
	}
	close V;	
}

sub readTags {
	my $line = shift;
	chomp $$line;
	my @row = split(/\t/, $$line);
	my @res = (@row[0 .. 6]);
	my @infoDat = split(/[;]/, $row[7]);
	my @formatDat1 = split(/:/, $row[8]);
	my %infoHsh = ();
	foreach my $i (@infoDat) {
		my @e = split(/=/, $i);
		$infoHsh{$e[0]} = defined($e[1])?$e[1]:"1";
	}
	
	foreach my $t (0 .. $#info) {
		push @res, defined($infoHsh{$info[$t]})?$infoHsh{$info[$t]}:".";
	}
	
	foreach my $s (9 .. $#row) {
		my %formatHsh = ();
		my @formatDat2 = split(/:/, $row[$s]);
		$formatHsh{$formatDat1[$_]} = $formatDat2[$_] foreach (0 .. $#formatDat1);

		foreach my $i (0 .. $#format) {
			push @res, defined($formatHsh{$format[$i]})?$formatHsh{$format[$i]}:".";	
		}
	}
	return \@res;
}
