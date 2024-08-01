#!/usr/bin/perl

# saf2H.pl
# usage: saf2H.pl <SAF text file of posterior probabilities>

use warnings;
use strict;

die(qq/
saf2H.pl <SAF text file of posterior probabilities>
\n/) if (!@ARGV && -t STDIN);

my $saf;
if (@ARGV) {
	open($saf, '<', $ARGV[0]) or die("Unable to open SAF file $ARGV[0]\n");
} else {
	die("No saf input detected\n") if (-t STDIN);
	$saf = \*STDIN;
}

my @sum = (0,0,0);
while(<$saf>) {
	chomp;
	my @tok = split(/\t/,$_);
	$sum[0] += exp($tok[2]);
	$sum[1] += exp($tok[3]);
	$sum[2] += exp($tok[4]);
}

my $total = $sum[0] + $sum[1] + $sum[2];
if ($total > 0) {
	my $h = $sum[1]/$total;
	print STDOUT "$h\n";
} else {
	print STDERR "Probabilities sum to zero\n";
}

exit;
