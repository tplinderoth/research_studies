#!/usr/bin/perl

# extractBiallelicSnps.pl

# notes:
# requires bcftools to be in user's path (no longer true as of v 2.0.0)
# Assumes input VCF indels are left-aligned and normalized: bcftools norm -f <fasta>

use warnings;
use strict;
use IO::Zlib;
use Getopt::Long;

my $version = '2.0.1';

my $minmaf = 0;
my $keep_overlap;
my $skip_af;

die(qq/
extractBiallelicSnps.pl v$version

extractBiallelicSnps.pl [parameters] <input vcf>
OR
cat VCF | extractBiallelicSnps.pl [parameters]

Input parameters:
--keep_overlap   Keep biallelic SNPs within deletions
--minmaf         FLOAT Only keep sites with minor allele frequency greater than FLOAT [$minmaf]
--skip_af        Do not check for INFO\/AF in range \(0,1\) for retaining SNPs

Notes:
* Assumes input VCF indels are left-aligned and normalized (bcftools norm -f <fasta>)
\n/) if (!@ARGV && -t STDIN);

# check whether bcftools executable is in user's path

#die("ERROR: Could not execute bcftools --> check that executable is in PATH") if system("bcftools > /dev/null 2>&1") < 0;
#chomp(my @bcftools_out = qx\bcftools version\);

#my $bcfv = $1 if ($bcftools_out[0] =~ /(\S+)$/);
#my $htsv = $1 if ($bcftools_out[1] =~ /(\S+)$/);

#my $depv = "bcftools-${bcfv}+htslib-${htsv}";

# parse user args

GetOptions('minmaf=f' => \$minmaf, 'keep_overlap' => \$keep_overlap, 'skip_af' => \$skip_af);
my %userargs = ('--minmaf' => $minmaf, '--keep_overlap' => $keep_overlap, '--skip_af' => $skip_af);

die("--minmaf should be a value in [0,0.5]\n") if ($minmaf < 0 || $minmaf > 0.5);

my $vcf = pop @ARGV if (@ARGV);

my $vcffh;
if ($vcf) {
	open($vcffh, '<', $vcf) or die("Unable to open VCF file $vcf: $!\n");
	sysread $vcffh, my $magic, 2;
	close $vcffh;

	if ($magic eq "\x1f\x8b") {
		# gzipped file
		$vcffh = new IO::Zlib;
		$vcffh->open($vcf, "rb");
	} else {
		open($vcffh, '<', $vcf);
	}
} else {
	die("No input VCF\n") if (-t STDIN);
	$vcffh = \*STDIN;
}

# declare variables to hold stats
my $nsites = 0; # number of sites in VCF
my $nsnps = 0; # total number SNPs
my $nsnp_sites; # total number sites containing SNPs
my $nsnp_widel = 0; # number of SNP sites within a deletion
my $nbisnps = 0; # number sites with biallelic snps
my $nbisnps_widel = 0; # number sites with biallelic snps within a deletion
my $ninsert = 0; # number of insertions
my $ndelete = 0; # number of deletions
my $nmask = 0; # number of times individual genotypes were masked
my $nout_sites = 0; # number SNPs in output VCF

# proces VCF header

my $datestr = localtime();
my $command = "##extractBiallelicSnpsCommand=<ID=extractBiallelicSnps.pl,Version=${version},Date=\"$datestr\",CommandLineOptions=\"";
foreach my $arg (keys %userargs) {
	$command .= "$arg $userargs{$arg} " if (defined $userargs{$arg});
}
$command .= "$vcf" if ($vcf);
$command =~ s/\s$//;
$command .= "\">\n";

# print pretty header

my @header = makeHeader($vcffh, $command);

foreach my $hline (@header) {
	print STDOUT $hline;
}

# Extract biallelic SNPs from VCF

$" = "\t";

my @indel = (0,0,''); # start and stop coordinates for indel
my ($end, $aa, $vartype, $chr);
while (my $line = <$vcffh>) {
	my @tok = split(/\s+/, $line);
	$nsites++;
	$chr = $tok[0];
	@indel = (0,0,'') if ($chr ne $indel[2]); # reset indel for new chromosome

	die("Multiple reference alleles found at $tok[0] $tok[1]") if ($tok[3] =~ /,/);
	die("Unexpecte allele in reference at $tok[0] $tok[1]") if ($tok[4] =~ /[^ACGT]i/);

	my $reflen = length($tok[3]);

	if ($reflen > 1) {
	# left-alignment means this is a deletion wrt reference
		# figure out if site is within a new deletion
		my $end = $tok[1] + ($reflen-1);
		if ($end > $indel[1]) {
			my $start = $tok[1]+1;
			$indel[0] = $tok[1]+1 if $start > $indel[1]; # this is the first base actually deleted
			$indel[1] = $end; # this is the last base actually deleted
			$indel[2] = $tok[0]; # chromosome
		}
	}
	my $widel = ($tok[1] >= $indel[0] && $tok[1] <= $indel[1] && $chr eq $indel[2]) ? 1 : 0;

	# classify the type of variant
	next if $tok[4] eq '.'; # skip nonvariable sites
	my @altarr = split(/,/, $tok[4]);
	my $issnp = 0;
	my $altn;
	my @snp_alleles = ();
	my $refbase = substr $tok[3], 0, 1;
	my $i = 1;

	foreach my $a (@altarr) {
		if ($a eq '*') {
			# do not record or keep star alleles
			$i++;
			next;
		}
		my $altlen = length($a);
		if ($altlen == $reflen) {
		# SNP
			$issnp = 1;
			my $altbase = substr $a, 0, 1;
			if ($altbase ne $refbase) {
				$altn = $i;
				push @snp_alleles, $altbase;
				$nsnps++;
			}
		} else {
		# indel
			$altlen > $reflen ? $ninsert++ : $ndelete++;
		}
		$i++;
	}

	$nsnp_sites++ if ($issnp);
	$nsnp_widel++ if ($issnp && $widel);

	# process biallelic SNP
	if (scalar(@snp_alleles) == 1) {
	# there is one alternate SNP allele (biallelic SNP site)
		$nbisnps++;

		# change representation of ref and kept alt allele in VCF line
		$tok[3] = $refbase;
		$altarr[$altn-1] = $snp_alleles[0];
		$tok[4] = join(',', @altarr);
		
		# mask individuals with nonref and alt snp alleles
		# could recalculate INFO subfields here while looping over genotypes but will leave that to bcftools for now

		# Update as of version 2.0.0: Shouldn't need to mask genotypes if using raw bcftools calls that have only been left-aligned
		# (no merging variants with norm -m). Multiallelic SNP entries are okay because these do not get printed.
		#if ($i > 2) {
		#	my $npl = ($i*($i+1)/2); # number of genotype likelihoods for nalleles alleles
		#	my $mask = "./.:" . "0," x ($i-1) . "0:0:0:" . "0," x ($npl-1) . "0"; #GT:AD:DP:GQ:PL, GQ=0 indicates masked individual
		#	
		#	for (my $j = 9; $j <= $#tok; $j++) {
		#		my ($a1, $a2) = ($1, $2) if ($tok[$j] =~ /^([.|\d]+)\/([.|\d]+):/);
                #
		#		if (defined $a1 && defined $a2) {
		#			if ($a1 eq '.' && $a2 eq '.') {
		#			next;
		#		}
		#		elsif (($a1 > 0 && $a1 != $altn) || ($a2 > 0 && $a2 != $altn)) {
		#			$tok[$j] = $mask;
		#			$nmask++;
		#		}
		#		} else {
		#			die("ERROR: Unrecognized genotype format at $tok[0] $tok[1]: $tok[$j]\n");
		#		}
		#
		#	}
		#}

		# adjust INFO ancestral alleles
		if ($tok[7] =~ /AA=([ACGTN*]+)/) {
			my $aa = $1;
			my $aa_base = substr $1, 0, 1;
			if ($aa_base eq $tok[3] || $aa_base eq $snp_alleles[0]) {
				$tok[7] =~ s/AA=$aa/AA=$aa_base/;
			} else {
				$tok[7] =~ s/AA=$aa/AA=N/
			}
		}

		# annotate INFO variant type
		$nbisnps_widel++ if ($widel);

		my $vt;
		if ($tok[7] =~ /VT=([^\s|;|=]+)/) {
			$vt = $1;
			my $vt_original = $vt;
			$vt .= ',snp' unless $vt =~ /snp/i;
			$vt .= ',widel' if ($widel && $vt !~ /widel/i);
			$tok[7] =~ s/VT=$vt_original/VT=$vt/;
		} else {
			$vt = 'snp';
			$vt .= ',widel' if ($widel);
			$tok[7] .= ";VT=$vt";
		}

		# print site
		unless ($widel && !$keep_overlap) {
			if (!$skip_af) {
				my $af;
				$af = $1 if ($tok[7] =~ /[;\s]+AF=([^;\s]+)/);
				if (!$af) {
					print STDERR "WARNING: No INFO/AF found for @tok[0..7]\n";
					next;
				}
				next unless ($af > 0.0 && $af < 1.0);
			}
			print STDOUT "@tok\n";
			$nout_sites++;
		}
	}

}

close $vcffh;

my $nvariants = $nsnps + $ndelete + $ninsert; # total number of variants

# print stats

print STDERR "Last contig processed: $chr\n";
print STDERR "Number VCF entries (lines without header): $nsites\n";
print STDERR "Number variants: $nvariants\n";
print STDERR "Number SNPs: $nsnps\n";
print STDERR "Number insertions: $ninsert\n";
print STDERR "Number deletions: $ndelete\n";
print STDERR "Number SNP sites: $nsnp_sites\n";
print STDERR "Number bialleleic SNP sites: $nbisnps\n";
print STDERR "Number SNP sites within a deletion: $nsnp_widel\n";
print STDERR "Number biallelic SNP sites within a deletion: $nbisnps_widel\n";
print STDERR "Number masked genotype entries: $nmask\n";
print STDERR "Number sites in biallelic SNP VCF: $nout_sites\n";

exit;

sub makeHeader {
	my ($fh, $command) = @_;

	# 1000 Genomes style VT annotation: ##INFO=<ID=VT,Number=.,Type=String,Description="indicates what type of variant the line represents">
	# widel denotes that the SNP is within a deletion
	my $vt_string = "##INFO=<ID=VT,Number=.,Type=String,Description=\"Denotes variant type\">\n";

	my @headorder = ('fileformat', 'reference', 'contig', 'INFO', 'FILTER', 'FORMAT', 'ALT', 'other'); # header order
	my %header = (fileformat => undef, ALT => undef, FILTER => undef, INFO => undef, FORMAT => undef, reference => undef, contig => undef, other => undef);
	my %seen = (vt => 0);
	my @headerlines;
	my $line;

	while (($line = readline($$fh)) =~ /^##/) {
		if ($line =~ /^##([^=]+)/i) {
			my $annotation = $1;

			if ($line =~ /INFO=<ID=VT,/) {
				# update VT annotation to ensure it is correct
				push @{$header{$annotation}}, $vt_string;
				$seen{vt} = 1;
			} elsif ($annotation =~ /^fileformat|ALT|FILTER|INFO|FORMAT|reference|contig/) {
				push @{$header{$annotation}}, $line;
			} else {
				push @{$header{other}}, $line;
			}
		}
	}

	push @{$header{INFO}}, "${vt_string}" if (!$seen{vt});
	push @{$header{other}}, $command if ($command);

	foreach my $tag (@headorder) {
		foreach (@{$header{$tag}}) {
			push @headerlines, $_;
		}
	}
	push @headerlines, $line;
	
	return @headerlines;
}
