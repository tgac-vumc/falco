#!/usr/bin/perl -w
#

my $mani = shift || "TruSeq_Amplicon_Cancer_Panel_Manifest_AFP1_PN15032433.txt";
open M, "<$mani";

my %probes = ();
while (<M>) {
	if (/\[Probes\]/) { readline(M); last; }
}

while (<M>) {
	if (/\[Targets\]/) {
		readline(M);
		last;
	}
	my @row = split(/\t/, $_);
	$probes{$row[2]} = [@row];
	#print STDERR join(":", @row[2, 9, 11]) . "\n"; sleep 1; 

}

while (<M>) {
	my @row = split(/\t/, $_);
	#print join("\t", @row[3, 4,5]) . "\n";
	my $ori = $probes{$row[0]}->[13];
	my $oriT = $row[6];
	my $s = $row[4] + length($probes{$row[0]}->[9]);
	my $e = $row[5] - length($probes{$row[1]}->[11]);
	if ($oriT eq "+") {
		$s = $row[4] + length($probes{$row[1]}->[11]);
		$e = $row[5] - length($probes{$row[0]}->[9]);
	}
	$row[4]--;
	$s--;
	print join("\t", @row[3, 4, 5], $row[0], 100, $oriT, $s, $e) . "\n";
}
