#!/usr/bin/perl

# saf2sfs.pl

use warnings;
use strict;

die(qq/
saf2sfs.pl <0=no fold | 1=fold> <SAF text file of posterior probabilities>
\n/) if ((!@ARGV && -t STDIN) || (!-t STDIN && !@ARGV) || (-t STDIN && @ARGV && scalar @ARGV < 2));

my $fold = $ARGV[0];
die("Invalid fold argument, must be 0|1\n") unless ($fold == 0 || $fold == 1);
my $saf;
if (@ARGV && scalar @ARGV == 2) {
	open($saf, '<', $ARGV[1]) or die("Unable to open SAF file $ARGV[0]\n");
} else {
	die("No saf input detected\n") if (-t STDIN);
	$saf = \*STDIN;
}

my @sum;
my @tok;
chomp(@tok = split(/\t/,<$saf>));
for (my $i = 2; $i <= $#tok; $i++) {
	push @sum, exp($tok[$i]);
}

while(<$saf>) {
	chomp;
	@tok = split(/\t/,$_);
	my $j = 0;
	for (my $i = 2; $i <= $#tok; $i++) {
		$sum[$j] += exp($tok[$i]);
		$j++;
	}
}

close $saf;

if ($fold) {
	my $n = (scalar @sum - 1)/2;
	for (my $i = 0; $i < $n; $i++) {
		my $k = $#sum - $i;
		$sum[$i] += $sum[$k];
		$sum[$k] = 0;
	}
}


print STDOUT "@sum\n";

exit;
