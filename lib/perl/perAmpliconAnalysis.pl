#!/usr/bin/perl -w

use strict;

my $bam = shift;
my $ref = shift;
my $manifest = shift;
my $base = shift || "output";
my $samtools = shift || `which samtools`;
my $fai = $ref . ".fai";
my $Qthreshold = 30;
my $QavgLim = 0;

open FAI, "<$fai";
my %idx = ();
print STDOUT "Reading fa index ... ";
while (<FAI>) {
        chomp;
        my @row = split(/\t/, $_);
        $idx{$row[0]} = [@row];
}
print STDOUT "done\n";

close FAI;
open REF, "<$ref";

my %stats = ();
my %manif = ();
my %targets = ();
my %probe = ();
my %hash = ();

open M, "<$manifest";
while (<M>) { last if (/\[Probes\]/) }
readline(M);
while (<M>) {
        chomp;
        if (/\[Targets\]/) {
                my $h = readline M;
                chomp $h;
                my @H = split(/\t/, $h);
                last;
        }
        my @row = split(/\t/, $_);
        $probe{$row[2]}{probe} = [@row];
        $probe{$row[2]}{ULSO} = $row[9];
        $probe{$row[2]}{DLSO} = $row[11];
}

while (<M>) {
        chomp;
        s/\s+$//;
        my @row = split(/\t/, $_);
        push @{$targets{$row[3]}{$_}}, $row[0] foreach ($row[4] .. $row[5]);
    #    $manif{$aTarget} = [@row];
	my $aTarget = $row[0] . "-" . $row[3] . ":". "$row[4]-$row[5]";
        if (not exists($stats{$aTarget}{target})) {
                $stats{$aTarget}{target} = [@row];
                $stats{$aTarget}{ln} = $row[5] - $row[4] + 1;
                $stats{$aTarget}{n} = 0;
                $stats{$aTarget}{sam} = [];
                $stats{$aTarget}{depth} = 0;
                $stats{$aTarget}{covs} = [];
                $stats{$aTarget}{chr} = $row[3];
                $stats{$aTarget}{start} = $row[4];
                $stats{$aTarget}{end} = $row[5];
		$stats{$aTarget}{ULSO} = $probe{$row[0]}{ULSO};
		$stats{$aTarget}{DLSO} = $probe{$row[0]}{DLSO};
                $stats{$aTarget}{assayStart} = $row[4] + length($probe{$row[0]}{ULSO});
                $stats{$aTarget}{assayEnd} = $row[5] - length($probe{$row[0]}{DLSO});
        }
#	print $row[0] . "\n";
	foreach my $i (0 .. 5) {
		foreach my $j (0 .. 5) {
#			print $row[4] + $i .":". ($row[5] - $row[4] + $j) .":". $row[0] . "\n";
#			print $row[4] + $i .":". ($row[5] - $row[4] - $j) .":". $row[0] . "\n";
#			print $row[4] - $i .":". ($row[5] - $row[4] + $j) .":". $row[0] . "\n";
#			print $row[4] - $i .":". ($row[5] - $row[4] - $j) .":". $row[0] . "\n";
			$hash{$row[3]}{$row[4] + $i}{$row[5] - $row[4] + $j} = $aTarget;
			$hash{$row[3]}{$row[4] + $i}{$row[5] - $row[4] - $j} = $aTarget;
			$hash{$row[3]}{$row[4] - $i}{$row[5] - $row[4] + $j} = $aTarget;
			$hash{$row[3]}{$row[4] - $i}{$row[5] - $row[4] - $j} = $aTarget;
		}
	}
#	sleep 1;
}
close M;

# samtools view -b 25MUT5.bam chr9:133738302-133738491 | samtools calmd - /gsdata/tmp/TSACP/ampliconPipeLine/ref/hg19/bwa-5/hg19.fa | head

#foreach my $amp (keys(%stats)) {

#	print $amp . "\n";
#	my @SAM = ();
	# Fetch reads
#	my $coord = $stats{$amp}{chr} .":".$stats{$amp}{start}."-".$stats{$amp}{end};
#	my $ln = $stats{$amp}{ln};
#	print $coord . " $ln\n";

# Open output 
open QC, ">$base.qc.ann.txt";
open QC2, ">$base.qc2.ann.txt";
open TARGET, ">$base.qc.targets.txt";

#TODO Remove stampy check
my $stampy = 0;
open SAMH, "$samtools view -H $bam |";
while (<SAMH>) {
	$stampy = 1 if (/PN:stampy/);
	print STDERR "STAMPY SPOTTED!!!" if (/PN:stampy/);
}
close SAMH;

open SAM, "$samtools view $bam |";
#open SAM, "$samtools view -b $bam | $samtools calmd - $ref |";
#	open SAM, "$samtools view -b $bam $coord  | $samtools calmd - $ref |";
while (<SAM>) {
	last unless (/^@/);
}
my $cnt = 0;
my %hsh = ();
my $pamp = ".";

#Print header

print QC "#" . join("\t", qw/chr pos amp dp ntRef ntVar nRef nVar nA nC nG nT nN ins del/) . "\n";
print QC2 "#" . join("\t", qw/read chr pos amp event/) . "\n";
print TARGET "#" . join("\t", qw/amp chr start end assayStart assayEnd depth/) . "\n";
while (<SAM>) {
	# Skip stampy reads? Nope, stampy is no longer needed
	#TODO Remove stampy check
#	next if ($stampy == 1 && ! /XP:Z:BWA/);
	chomp;
	my @row = split(/\t/, $_);
	my $l = 0;
	my $cig = $row[5];
	$cig =~ s/(\d+)[MD]/$l+=$1/eg;
		
	my $amp = $hash{$row[2]}{$row[3]}{$l} || next;
	#print "$row[2] $row[3] $l $amp\n"; 

	push @{$stats{$amp}{sam}}, [@row];
	$stats{$amp}{depth}++;
	if (($pamp ne ".") && ($amp ne $pamp)) {
		print STDOUT "Flushing $pamp $stats{$pamp}{depth}\n";
#		print TARGET join("\t", $pamp, map { $stats{$pamp}{$_} } qw/chr start end assayStart assayEnd depth/) . "\n";
		my $data = &call($stats{$pamp}{sam});
		&flush($data, $stats{$pamp}{chr}, $pamp, $stats{$pamp}{assayStart}, $stats{$pamp}{assayEnd});
		$stats{$pamp}{sam} = [];
	}
	$pamp = $amp;
	$cnt++;
}

print STDOUT "Flushing $pamp $stats{$pamp}{depth}\n";
my $data = &call($stats{$pamp}{sam});
&flush($data, $stats{$pamp}{chr}, $pamp, $stats{$pamp}{assayStart}, $stats{$pamp}{assayEnd});
$stats{$pamp}{sam} = [];

close QC;
close QC2;
foreach my $pamp (keys %stats) {
		print TARGET join("\t", $pamp, map { $stats{$pamp}{$_} || "NA"} qw/chr start end assayStart assayEnd depth/) . "\n";
}
close TARGET;
	

#	print join("\n", map {$hsh{$_} . ":" . $_} keys(%hsh)) . "\n";
#	print $cnt . "\n"; sleep 1;
	close SAM;
#	$stats{$amp}{depth} = $cnt;
#	my $data = &call(\@SAM);
#	&flush($data, $stats{$amp}{chr}, $amp);

#}

#foreach my $amp (keys(%stats)) {

#	print $amp . ":" . scalar(@{$stats{$amp}{sam}}) . "\n";
#}

sub call {
	my $SAM = shift;

	my $ppos = 0;
	my $pchr = "";
	my %data = ();
	my $refSeq = "";
	my $c = 0;
	my $col_chr = 2;
	my $col_read = 9;
	my $col_qual = 10;
	my $col_pos = 3;
	my $col_cigar = 5;
	my %cooc = ();
	my $cnt = 1;
	foreach my $r (0 .. $#$SAM) {
		my @row = @{$SAM->[$r]};
		my $read = $row[$col_read];
		my $qual = $row[$col_qual];
		my $pos = $row[$col_pos];
		my @cig = ();
		my $aln = 0;

		# Calculate mean read quality and skip if it's below threshold
		if ($QavgLim > 0) {
			my $QSum = 0;
			my $Qlen = length($qual);

			$QSum += ord($qual) foreach 1 .. $Qlen;

			my $Qavg = $QSum / $Qlen - 33;
			next if ($Qavg < $QavgLim);
		}

		while ($row[$col_cigar] =~ /(\d+)(\D)/g) {
			push @cig, [$1, $2];
			$aln += $1;
		}
		if ($ppos != $pos) {
			my $offset = $idx{$row[$col_chr]}->[2]; # Byte offset for chromosome
			my $nl = int (($row[$col_pos] - 1) / $idx{$row[$col_chr]}->[3]); # Newline bytes to offset
			$offset += $row[$col_pos] + $nl - 50 - 1; # Read in 50 nt prior to start of amplicon
			seek REF, $offset, 0;
			read REF, $refSeq, 350;
			$refSeq =~ s/\s//g;
			$data{$row[$col_chr]}{$row[$col_pos] + $_ - 50 + 1}{"refNT"} = substr($refSeq, $_, 1) foreach (0 .. 300);
			$ppos = $pos;
		}

		foreach my $cigE (@cig) {
			my $n = $cigE->[0];
			if ($cigE->[1] eq "I") {
				my $insert = $data{$row[$col_chr]}{$pos - 1}{"refNT"};
				$insert .= substr($read, 0, $n, "");
				my $qinsert .= substr($qual, 0, $n, "");
				push @{$data{$row[$col_chr]}{$pos - 1}{"IR"}{$insert}}, $r;
				$data{$row[$col_chr]}{$pos - 1}{"I"}{$insert}++;
			}
			elsif ($cigE->[1] eq "M") {
				my $mseq = substr($read, 0, $n, "");
				my $mqual = substr($qual, 0, $n, "");
				foreach my $nt (0 .. (length($mseq) - 1)) {
					$data{$row[$col_chr]}{$pos + $nt}{"DP"}++;
					my $ntq = substr($mqual, $nt, 1);
					my $ntn = substr($mseq, $nt, 1);
					if (ord($ntq) - 33  < $Qthreshold) {
						$ntn = "N";
					}
					$data{$row[$col_chr]}{$pos + $nt}{"NT"}{$ntn}++;
					push @{$data{$row[$col_chr]}{$pos + $nt}{"NTR"}{$ntn}}, $r if (uc($ntn) ne uc($data{$row[$col_chr]}{$pos + $nt}{refNT}));
				}
				$pos += $n;

			}
			elsif ($cigE->[1] eq "D") {
				my $deletion = join("", map { $data{$row[$col_chr]}{$pos + $_ -1}{"refNT"} } (0 .. ($n)));
				push @{$data{$row[$col_chr]}{$pos - 1}{"DR"}{$deletion}}, $r;
				$data{$row[$col_chr]}{$pos - 1}{"D"}{$deletion}++;
				$pos += $n;
			}
		}
		$cnt++;
	}
	#my $thresh = scalar(@$SAM) * .05;
	#&coocCalc(\%cooc, \%data, $thresh);	
	return \%data;
}

sub coocCalc {
	my $cooc = shift;
	my $data = shift;
	my $thresh = shift;
	my %rev = ();
	# Filter low frequency variants 
	foreach my $k (keys(%$cooc)) {
		if (scalar(keys(%{$cooc->{$k}})) < $thresh) {
				delete($cooc->{$k});
				next;
		}
		print STDERR "$k : $thresh : " . join(":", sort {$a cmp $b} keys(%{$cooc->{$k}})) . "\n";
	}
	my @muts = keys(%{$cooc});
	foreach my $i (0 .. $#muts) {
#		my %hsh = map { $_ => 0 } keys(%{$cooc{$muts[$i]}});
		foreach my $j (($i + 1) .. $#muts) {
			my $score = 0;
			foreach my $m (keys(%{$cooc->{$muts[$i]}})) {
				if (exists($cooc->{$muts[$j]}{$m})) {
					$score++;	
				}
			}
			print STDERR $muts[$i] .":". $muts[$j] . ":" . $score . "\n"; sleep 1;
		}	
	}		

	#{read}{mutation}	
	#{mutation}{read}	
	sleep 1;
}

sub flush {
        my $data = shift;
        my $chr = shift;
	my $amp = shift;
	my $s = shift;
	my $e = shift;
        return unless (exists($data->{$chr}));
        print STDOUT "printing ... ";
	my %reads = (); # read->pos->nt
#        foreach my $p (sort {$a <=> $b} keys(%{$data->{$chr}})) {
#        	unless (($p >=  $s) && ($p <= $e)) {
#			next;	
#		}
	foreach my $p ($s .. $e) {
	        my $posdat = $data->{$chr}{$p} || next;
	        next unless (exists($posdat->{DP}));
                my $refNT = uc($posdat->{refNT});
                my @varOrd = sort { $posdat->{NT}{$b} <=> $posdat->{NT}{$a} } grep {$_ ne "N"} grep { $_ ne $refNT  } keys(%{$posdat->{NT}});
                my $seqNT = join("/", @varOrd);
                #my $refCNT = scalar(@{$posdat->{NT}{$refNT} || []});
		my $refCNT = $posdat->{NT}{$refNT} || 0;
		my $varCNT = 0;
		#$varCNT += scalar(@{$posdat->{NT}{$_}}) foreach (@varOrd);
		#$varCNT += $posdat->{NT}{$_} foreach (@varOrd);
		$varCNT = ($varOrd[0])?$posdat->{NT}{$varOrd[0]}:0;
		my $Istr = ".";
		my $Dstr = ".";
		if ($varOrd[0]) {
			$reads{$_}{$p} = "$refNT/$varOrd[0]" foreach (@{$posdat->{NTR}{$varOrd[0]}});
		}
		if (exists($posdat->{I})) {
			$Istr = join("|", map { $_ . ":" . $posdat->{I}{$_} } keys(%{$posdat->{I}}));
			foreach my $I (keys(%{$posdat->{IR}})) {
				$reads{$_}{$p} = "$refNT/$I" foreach (@{$posdat->{IR}{$I}});
			}
		}
		if (exists($posdat->{D})) {
			$Dstr = join("|", map { $_ . ":" . $posdat->{D}{$_} } keys(%{$posdat->{D}}));
			foreach my $D (keys(%{$posdat->{DR}})) {
				$reads{$_}{$p} = "$D/$refNT" foreach (@{$posdat->{DR}{$D}});
			}
		}
		#$Dstr = join("|", map { $_ . ":" . $posdat->{D}{$_} } keys(%{$posdat->{D}})) if exists($posdat->{D});
		print QC join("\t", $chr, $p, $amp,
				$posdat->{DP} || 0,
				$posdat->{refNT} || ".",
				$seqNT || ".",
				$refCNT, $varCNT,
				#(map { scalar(@{$posdat->{NT}{$_} || []})} qw(A C G T N)), $Istr, $Dstr,
				(map { $posdat->{NT}{$_} || 0 } qw(A C G T N)), $Istr, $Dstr,
			  ) . "\n";
		# delete $data->{$chr}{$p};
	}
	delete $data->{$chr};
	foreach my $r (sort { $a <=> $b } keys(%reads)) {
		print QC2 join("\t", $r, $chr, $_, $amp, $reads{$r}{$_}) . "\n" foreach (sort {$a <=> $b} keys(%{$reads{$r}})); 
	}
	print STDOUT "done\n";
}

