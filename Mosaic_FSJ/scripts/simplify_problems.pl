#!/usr/bin/perl

# simplify_problems.pl <input file>

use warnings;
use strict;

die(qq/
simplify_problems.pl <input file>
\n/) unless @ARGV && scalar @ARGV == 1;

open(my $fh, '<', $ARGV[0]) or die("Unable to open file $ARGV[0]: $!\n");

print STDOUT "SOURCE file: 'JAY ID DATA WITH PARENTS through 2021.xlsm', tab: 'JAYID'\n";
print STDOUT "METADATA file: 'M4_FSJ_Blood_Samples.xlsx', tab: 'SURVIVORSHIP_WNESTING'\n\n";

while(<$fh>) {
	chomp;
	if ($_ eq "### BOTH DAM AND SIRE ###") {
		print STDOUT "$_\n";
		print STDOUT "ID\tSOURCE_DAM_OF\tSOURCE_SIRE_OF\tMETA_DAM_OF\tMETA_SIRE_OF\n";
		while ((my $line = <$fh>) =~ /^\S+/) {
			chomp($line);
			#if ($line =~ /\[SOURCE/ && $line =~ /\[SURVIVORSHIP_WNESTING/) {
			my $id = $1 if $line =~ /^(\S+)/;
			my $sd = $1 if $line =~ /\[SOURCE:\sdam_of:\s([^\]]+)\ssire_of/;
			$sd = "NA" unless $sd && $sd =~ /\S/;
			$sd =~ s/\s/,/g;
			my $ss = $1 if $line =~ /\[SOURCE:[^\]]+sire_of:\s([^\]]+)\]/;
			$ss = "NA" unless $ss && $ss =~ /\S/;
                        $ss =~ s/\s/,/g;
			my $md = $1 if $line =~ /\[SURVIVORSHIP_WNESTING:\sdam_of:\s([^\]]+)\ssire_of/;
			$md = "NA" unless $md && $md =~ /\S/;
                        $md =~ s/\s/,/g; 
			my $ms = $1 if $line =~ /\[SURVIVORSHIP_WNESTING:[^\]]+sire_of:([^\]]+)\]/;
                        $ms = "NA" unless $ms && $ms =~ /\S/;
			$ms =~ s/^\s+|\s+$//g;
                        $ms =~ s/\s/,/g;
			print STDOUT "$id\t$sd\t$ss\t$md\t$ms\n";
			#}
			last if eof($fh);
		}
		print STDOUT "\n";
	} elsif ($_ eq '### DIFFERENT PARENTS ###') {
		print STDOUT "$_\n";
		print STDOUT "ID\tSOURCE_MOM\tSOURCE_DAD\tMETA_MOM\tMETA_DAD\n";
		while ((my $line = <$fh>) =~ /^\S+/) {
			chomp($line);
			if ($line =~ /\[SOURCE/ && $line =~ /\[SURVIVORSHIP_WNESTING/) {
				my $id = $1 if $line =~ /^(\S+)/;
				my $smom = $1 if $line =~ /\[SOURCE:\smother:\s([\S]+)\sfather/;
				$smom = "NA" unless $smom && $smom =~ /\S/ && lc($smom) ne "unknown";
				my $sdad = $1 if $line =~ /\[SOURCE:[^\]]+father:\s([^\]]+)\]/;
				$sdad = "NA" unless $sdad && $sdad =~ /\S/ && lc($sdad) ne "unknown";
				my $mmom = $1 if $line =~ /\[SURVIVORSHIP_WNESTING:\smother:\s([\S]+)\sfather/;
				$mmom = "NA" unless $mmom && $mmom =~ /\S/ && lc($mmom) ne "unknown";
				my $mdad = $1 if $line =~ /\[SURVIVORSHIP_WNESTING:[^]]+father:\s([^\]]+)\]/;
				$mdad = "NA" unless $mdad && $mdad =~ /\S/ && lc($mdad) ne "unknown";
				print STDOUT "$id\t$smom\t$sdad\t$mmom\t$mdad\n" if ($smom ne $mmom || $sdad ne $mdad);
			}
			last if eof($fh);
		}
		print STDOUT "\n";
	} elsif ($_ eq '### DIFFERENT SEX ###') {
		print STDOUT "$_\n";
		print STDOUT "ID\tSOURCE_SEX\tMETA_SEX\n";
		while ((my $line = <$fh>) =~ /^\S+/) {
			chomp($line);
			if ($line =~ /\[SOURCE/ && $line =~ /\[SURVIVORSHIP_WNESTING/) {
				chomp($line);
				my $id = $1 if $line =~ /^(\S+)/;
				my $ssex = $1 if $line =~ /\[SOURCE\ssex:\s([^\]]+)/;
				my $msex = $1 if $line =~ /\[SURVIVORSHIP_WNESTING\ssex:\s([^\]]+)/;
				print "$id\t$ssex\t$msex\n";
			}
			last if eof($fh);
		}
		print STDOUT "\n";
	} elsif ($_ eq "### MISSING FROM SOURCE \(file: 'JAY ID DATA WITH PARENTS through 2021.xlsm', tab: 'JAYID'\) ###") {
		print STDOUT "$_\n";
		while ((my $line = <$fh>) =~ /^\S+/) {
			print STDOUT $line;
			last if eof($fh);
		}
		print STDOUT "\n";
	} elsif ($_ eq "### MISSING FROM SURVIVORSHIP_WNESTING METADATA \(file: 'M4_FSJ_Blood_Samples.xlsx', tab: 'SURVIVORSHIP_WNESTING'\) ###") {
		print STDOUT "$_\n";
		while ((my $line = <$fh>) =~ /^\S+/) {
			print STDOUT $line;
			last if eof($fh);
		}
	}
}

close $fh;

exit;
