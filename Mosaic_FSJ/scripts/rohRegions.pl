#!/usr/bin/perl

# rohRegions.pl

use warnings;
use strict;

my $minq = 20;

die(qq/
rohRregions.pl [minimum site quality for FROH estimation, default:$minq] <bcftools per site ROH file> <outfile name>
\n/) if (!@ARGV || scalar @ARGV < 2);

my $outname = pop @ARGV;
my $inname = pop @ARGV;
$minq = $ARGV[0] if @ARGV;

open(my $ifh, '<', $inname) or die("Unable to open per site ROH file $inname: $!\n");
open(my $ofh, '>', $outname) or die("Unable to open output regions file $outname: $!\n");

my ($nsites, $nsites_hiq, $nauto, $nauto_hiq, $prevstate, $rohlen, $avgq, $prevscaff) = (0, 0, 0, 0, 0, 0, 0, '');
my @reg = ('', 0, 0, 0, 0); # scaffold, start, end, quality, number SNPs in ROH region
$" = "\t";

print $ofh "SCAFFOLD\tROH_START\tROH_END\tLENGTH\tAVERAGE_QUALITY\n";

while (<$ifh>) {
	next if $_ !~ /^ROH/;
	chomp;
	my @tok = split(/\t/, $_);
	
	if ($tok[2] ne $prevscaff) {
		# new scaffold
		if ($prevstate == 1) {
			# write ROH region if last scaffold ended in a ROH
			$rohlen = $reg[2] - $reg[1] + 1;
			$avgq = sprintf("%.6f",$reg[3]/$reg[4]);
			print $ofh "@reg[0..2]\t$rohlen\t$avgq\n";
			# reset region array
			@reg = ('', 0, 0, 0, 0);
			# clear previous state
			$prevstate = 0;
		}
	}
	
	if ($tok[4] == 1) {
		# new site is autozygous
		$nauto++;
		$nauto_hiq++ if $tok[5] >= $minq;
		if ($prevstate == 0) {
			# new ROH region
			$reg[0] = $tok[2];
			$reg[1] = $tok[3];
			$reg[2] = $tok[3];
			$reg[3] = $tok[5];
			$reg[4] = 1;
		} else {
			$reg[2] = $tok[3]; # update end position
			$reg[3] += $tok[5];
			$reg[4]++;
		}
	} else {
		# new site is not autozygous
		if ($prevstate == 1) {
			# end of ROH region
			$rohlen = $reg[2] - $reg[1] + 1;
			$avgq = sprintf("%.6f",$reg[3]/$reg[4]);
			print $ofh "@reg[0..2]\t$rohlen\t$avgq\n";
			# reset region array
			@reg = ('', 0, 0, 0, 0);
		}
	}
	$nsites++;
	$nsites_hiq++ if $tok[5] >= $minq;
	$prevstate = $tok[4];
	$prevscaff = $tok[2];
}

if ($prevstate == 1) {
	$rohlen = $reg[2] - $reg[1] + 1;
	$avgq = sprintf("%.6f",$reg[3]/$reg[4]);
	print $ofh "@reg[0..2]\t$rohlen\t$avgq\n";
}

close $ifh;
close $ofh;

my $pauto = $nsites > 0 ? $nauto/$nsites : 'NaN';
$pauto = sprintf("%.8f", $pauto); # format precision
my $pauto_hiq = $nsites_hiq > 0 ? $nauto_hiq/$nsites_hiq : 'NaN';
$pauto_hiq = sprintf("%.8f", $pauto_hiq); # format precision
print STDOUT "FROH\tFROH_HIQUAL\tN_SITES\tN_AUTOZYGOUS\tN_SITES_HIQUAL\tN_AUTOZYGOUS_HIQUAL\n$pauto\t$pauto_hiq\t$nsites\t$nauto\t$nsites_hiq\t$nauto_hiq\n";

exit;
