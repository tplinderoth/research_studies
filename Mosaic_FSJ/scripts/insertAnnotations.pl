#!/usr/bin/perl

# insertAnnotations.pl

use warnings;
use strict;
use Getopt::Long;
use IO::Zlib;
use List::MoreUtils qw(uniq);
#use warnings FATAL => 'all';

my $version = '1.1.0'; # updated to work with FSJ mosaic data
my $alfile = undef;
my $alfields = undef;
my $anctype = 'parsimony';
my $dpbounds = undef;
my $hetbound = undef;
my $vcf = undef;
my $rmfilter = undef;
my $fstart = 4;
my $bed = undef;
my $bedlist = undef;
my $genorep = undef;
my $overwrite;

die(qq/
insertAnnotations.pl version $version

insertAnnotations.pl [options] <vcf file>
OR
cat <vcf file> | insertAnnotations.pl [options]

--alfile    File of outgroup alleles in pafAlleles format
--alfields  Column numbers of alleles to use from alleles file, 1-indexed and ','-separated [$fstart...]
--anctype   Infer ancestral allele using 'parsimony' or as the 'major' allele among outgroups [$anctype]
--dpbounds  LowDP and HighDP bounds to annotate FILTER field, format '<INT lower bound>,<INT upper bound>'
--hetbound  Lower p-value bound for excess heterozygosity test \(uses INFO\/ExcHet field\)
--rmfilter  ','-delimited list of FILTER annotations to remove
--bed       Bed format file with FILTER annotations to add in the 4th column
--bedlist   ','-delimited list of FILTER annotations to ignore from bed
--overwrite Replace all existing filters, otherwise new filters are appended
--genorep   File of sample subset indexes for appending genotype representation to the INFO field

Notes:
*Assumes the same contig order in VCF and allele files.

*For inferring ancestral alleles with parsimony it is assumed that outgroups in the pafAlleles file 
 are ordered with increasing divergence from the ingroup (greater column number = greater divergence).
 Alternatively, the outgroup columns can be passed in order of their increasing divergence from the
 ingroup using the --alfields option.

*Parsimony currently only implemented for 3 outgroup case.

*Particular FILTER annotations will be updated if the function that applies them (--dpbounds, --bed)
 is called.

*Bed must be sorted in same order as VCF.

*--genorep file must have columns <group name> <comma-delimited, 0-based indexes of samples in VCF> <min GQ> <min DP>.
 Lines with '#' are ignored (so preceed a header with # for example).
\n/) if (-t STDIN && !@ARGV);

GetOptions('alfile=s' => \$alfile, 'alfields=s' => \$alfields, 'dpbounds=s' => \$dpbounds, 'hetbound=f' => \$hetbound, 'rmfilter=s' => \$rmfilter, 
'anctype=s' => \$anctype, 'bed=s' => \$bed, 'bedlist=s' => \$bedlist, 'overwrite' => \$overwrite, 'genorep=s' => \$genorep);

# this is for printing the command to the VCF header
my %userargs = ('--alfile' => $alfile, '--alfields' => $alfields, '--dpbounds' => $dpbounds, '--hetbound' => $hetbound, '--rmfilter' => $rmfilter, 
'--anctype' => $anctype, '--bed' => $bed, '--bedlist' => $bedlist, '--overwrite' => $overwrite, '--genorep' => $genorep);
# disable printing commands if not applicable
$userargs{"--anctype"} = undef if (!$alfile);
$userargs{"--bedlist"} = undef if (!$bedlist);
$userargs{"--genorep"} = undef if (!$genorep);

# read in VCF

$vcf = pop @ARGV;
my $vcffh;
if ($vcf) {
	open($vcffh, '<', $vcf) or die("Unable to open VCF $vcf: $!\n");
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


# initalize ancestral allele inference

my ($alfh, $nalleles, $aamethod, @outal, @col, @parstable, %tree_counts);
my ($l1, $l2, $chr, $nsites) = (undef, undef, "", 0);

if ($alfile) {
	open($alfh, '<', $alfile) or die("Could not open allele file $alfile: $!\n");

	if ($alfields) {
		@col = map {die("Argument to --alfields must be ','-delimited integers > 0") if $_ =~ /\D/; $_ - 1} split(/,/, $alfields);
                $nalleles = scalar @col;
	} else {
		$nalleles = split(/\s+/, <$alfh>)-3;
		@col = ($fstart-1) .. ($fstart + $nalleles - 2);
	}

        $l1 = <$alfh>; # first site of new chr
        die("No sites in allele file $alfile\n") if (!$l1);

	if (lc($anctype) eq 'parsimony') {
		$aamethod = 1;
		@parstable = makeLookup(\%tree_counts);
		treeCounter(\%tree_counts, $nalleles, 0);
	} elsif (lc($anctype) eq 'major') {
		$aamethod = 2;
	} else {
		die("Invalid --anctype argument: $anctype\n");
	}
} else {
	$anctype = undef;
}

# initialize genotype representation annotations

my (%grp, %grpseen, @grpid, %grpcounts);

open(my $grpfh, '<', $genorep) or die("Unable to open sample group index file $genorep: $!\n");
while(<$grpfh>) {
	if ($_ !~ /^#/) {
		chomp;
		my @tok = split(/\s+/, $_);
		if (! exists $grpseen{$tok[0]}) {
			push @grpid, $tok[0];
			$grpseen{$tok[0]} = 1;
			$grpcounts{$tok[0]}{gq} = 0;
			$grpcounts{$tok[0]}{dp} = 0;
		}
		die("Group info file must have four columns <group ID> <VCF sample indexes> <min GQ> <min DP>: $!\n") if scalar @tok != 4;
		foreach my $idx (split(/,/,$tok[1])) {
			my @idx2 = $idx =~ /(\d+)-(\d+)/ ? ($1..$2) : $idx;
			foreach my $i (@idx2) {
				$grp{$i}{id} = $tok[0];
				$grp{$i}{gq} = $tok[2];
				$grp{$i}{dp} = $tok[3];
			}
		}
	}
}
close $grpfh;

# initialize filter annotation

my (@dpcutoff, %delfilter, $bedfh);

if ($dpbounds) {	
	@dpcutoff = split(/,/, $dpbounds);
	if (scalar(@dpcutoff) != 2) {
		die("ERROR: --dpbounds requires 2 values, a lower and an upper coverage threshold\n");
	} elsif ($dpcutoff[0] >= $dpcutoff[1]) {
		die("ERROR: --dpbounds upper bound must be greater than the lower bound value\n");
	}
	
	$delfilter{HighDP} = 1;
	$delfilter{LowDP} = 1;
}

if ($hetbound) {
	die ("--hetbound outside [0,1] interval\n") if ($hetbound < 0 || $hetbound > 1);
	$delfilter{HighHet} = 1;
}

my %bedfilter = (LowCov => 1, ExcessCov => 1, LowMQ => 1, ExcessMQ0 => 1, NoCov => 1, LowBQ => 1, ExcessBQ0 => 1);
if ($bed) {
	if ($bedlist) {
		$bedlist .= "," if ($bedlist !~ /,$/);
		$bedlist =~ s/,/\|/g;
	}
	open($bedfh, '<', $bed) or die("Unable to open bed file $bed: $!\n");
	foreach (keys %bedfilter) {
		if ($bedlist && $bedlist =~ /$_\|/) {
			delete $bedfilter{$_};
			next;
		}
		$delfilter{$_} = 1;
	}
}

if ($rmfilter) {
	foreach my $filter_anno (split(/,/,$rmfilter)) {
		$delfilter{$filter_anno} = 1;
	}
}

# print pretty VCF header

my $datestr = localtime();
my $command = "##insertAnnotationsCommand=<ID=insertAnnotations.pl,Version=${version},Date=\"$datestr\",CommandLineOptions=\"";
foreach my $arg (keys %userargs) {
	$command .= "$arg $userargs{$arg} " if (defined $userargs{$arg});
}
$command .= "$vcf" if ($vcf);
$command =~ s/\s$//;
$command .= "\">\n";

# Should make subsequent versions take a header annotation file to generalize this for different
# filtering parameters.

my ($grp_gq_string, $grp_dp_string);

my $aa_string = "##INFO=<ID=AA,Number=1,Type=String,Description=\"Ancestral Allele\">\n";
if ($genorep) {
	$grp_gq_string = "##INFO=<ID=REPGQ,Number=" . scalar(@grpid) . ",Type=Integer,Description=\"Number of genotypes with GQ of 15 or above for classes " . join(', ',@grpid) . "\">\n";
	$grp_dp_string = "##INFO=<ID=REPDP,Number=" . scalar(@grpid) . ",Type=Integer,Description=\"Number of genotyped samples with DP of 3 or greater for classes " . join(', ',@grpid) . "\">\n";
}
my $lowdp_string = "##FILTER=<ID=LowDP,Description=\"Site DP less than genome-wide autosomal or sex chromosome median DP -25%\">\n";
my $highdp_string = "##FILTER=<ID=HighDP,Description=\"Site DP greater than genome-wide autosomal or sex chromosome median DP +25%\">\n";
my $highhet_string = "##FILTER=<ID=HighHet,Description=\"ExcHet below genome-wide 0.02 quantile\">\n";
my $lowcov_string = "##FILTER=<ID=LowCov,Description=\"Total site depth below genome-wide autosomal or sex chromomsome 0.05 quantile\">\n";
my $excesscov_string = "##FILTER=<ID=ExcessCov,Description=\"Total site depth above genome-wide autosomal or sex chromosome 0.95 quantile\">\n";
my $lowmq_string = "##FILTER=<ID=LowMQ,Description=\"Root mean square map quality below 35\">\n";
my $mqzero_string = "##FILTER=<ID=ExcessMQ0,Description=\"Fraction of reads with map quality zero exceeds 2x the genome-wide average\">\n";
my $lowbq_string = "##FILTER=<ID=LowBQ,Description=\"Root mean square base quality below 20\">\n";
my $bqzero_string = "##FILTER=<ID=ExcessBQ0,Description=\"Fraction of reads with base quality zero exceeds 2x the genome-wide average\">\n";
my $nocov_string = "##FILTER=<ID=NoCov,Description=\"No mapped reads\">\n";

my @headorder = ('fileformat', 'reference', 'contig', 'INFO', 'FILTER', 'FORMAT', 'ALT', 'other'); # header order

my %header = (fileformat => undef, ALT => undef, FILTER => undef, INFO => undef, FORMAT => undef, reference => undef, contig => undef, other => undef);

my %seen = (aa => 0, lowdp => 0, highdp => 0, highhet => 0);

my $hline;
my $filter;
while (($hline = <$vcffh>) =~ /^##/) {
	if ($hline =~ /^##([^=]+)/i) {
		my $annotation = $1;

		if (($rmfilter || $bed) && $hline =~ /^##FILTER=<ID=([^,]+)/) {
			next if exists $delfilter{$1}
		}

		if ($alfile && $hline =~ /INFO=<ID=AA/) {
			# update AA annotation to ensure it is correct
			push @{$header{$annotation}}, $aa_string;
			++$seen{aa};
		} elsif ($genorep && $hline =~ /INFO=<ID=REPGQ>/) {
			# update REPGQ annotation
			push @{$header{$annotation}}, $grp_gq_string;
			++$seen{repgq};
		} elsif ($genorep && $hline =~ /INFO=<ID=REPDP>/) {
			# update REPDP annotation
			push @{$header{$annotation}}, $grp_dp_string;
			++$seen{repdp};
		} elsif ($dpbounds && $hline =~ /FILTER=<ID=LowDP/) {
			# update LowDP annotation
			push @{$header{$annotation}}, $lowdp_string;
			++$seen{lowdp};
		} elsif ($dpbounds && $hline =~ /FILTER=<ID=HighDP/) {
			# update HighDP annotation
			push @{$header{$annotation}}, $highdp_string;
			++$seen{highdp};
		} elsif ($hetbound && $hline =~ /FILTER=<ID=HighHet>/) {
			# update HighHet annotation
			push @{$header{$annotation}}, $highhet_string;
			++$seen{highhet};
		} elsif ($annotation =~ /^fileformat|ALT|FILTER|INFO|FORMAT|reference|contig/) {
			push @{$header{$annotation}}, $hline;
		} else {
			push @{$header{other}}, $hline;
		}
	}
}

push @{$header{INFO}}, $aa_string if ($alfile && !$seen{aa});
push @{$header{INFO}}, $grp_gq_string if ($genorep && !$seen{repgq});
push @{$header{INFO}}, $grp_dp_string if ($genorep && !$seen{repdp});
push @{$header{FILTER}}, $lowdp_string if ($dpbounds && !$seen{lowdp});
push @{$header{FILTER}}, $highdp_string if ($dpbounds && !$seen{highdp});
push @{$header{FILTER}}, $highhet_string if ($hetbound && !$seen{highhet});
push @{$header{FILTER}}, $nocov_string if ($bedfilter{NoCov});
push @{$header{FILTER}}, $lowcov_string if ($bedfilter{LowCov});
push @{$header{FILTER}}, $excesscov_string if ($bedfilter{ExcessCov});
push @{$header{FILTER}}, $lowmq_string if ($bedfilter{LowMQ});
push @{$header{FILTER}}, $mqzero_string if ($bedfilter{ExcessMQ0});
push @{$header{FILTER}}, $lowbq_string if ($bedfilter{LowBQ});
push @{$header{FILTER}}, $bqzero_string if ($bedfilter{ExcessBQ0});
push @{$header{other}}, $command;


# make matrix for quick genorep logging
my @grpmat;
if ($genorep) {
	my $fieldidx = 0;
	foreach (split(/\s+/, $hline)) {
		push @grpmat, [0,0,0,0];
		if ($fieldidx > 8) {
			my $sampidx = $fieldidx - 9;
			if (exists $grp{$sampidx}) {
				$grpmat[$fieldidx][0] = 1;
				$grpmat[$fieldidx][1] = $grp{$sampidx}{id};
				$grpmat[$fieldidx][2] = $grp{$sampidx}{gq};
				$grpmat[$fieldidx][3] = $grp{$sampidx}{dp};
			}
		}
		$fieldidx++;
	}
}

# print header info

foreach my $tag (@headorder) {
	foreach (@{$header{$tag}}) {
		print STDOUT $_;
	}
}

print STDOUT $hline;

# process VCF sites

my %filtercounts;
map { $filtercounts{$1} = 0 if $_ =~ /ID=([^,]+)/ } @{$header{FILTER}};
my %filterconfigs;

my $vcfchr;
my @bedtok = ('',0,0);
my $site_count = 0;
my $fail_counts = 0;
my $aa_count = 0;
my @tok;
my %missing_chr;

$" = "\t";

while (<$vcffh>) {
	@tok = split(/\s+/, $_);
	$site_count++;

	## update bedgraph region - assumes VCF and bed are sorted in the same order
	if ($bed) {
		$vcfchr = uc($tok[0]);
		# update chromosome
		while (!eof($bedfh) && $bedtok[0] ne $vcfchr) { bedsplit(\@bedtok, $bedfh); }
		if (eof($bedfh) && $bedtok[0] ne $vcfchr) {
			if (!exists $missing_chr{$tok[0]}) {
				print STDERR "WARNING: Reached the end of bedgraph file before finding $tok[0]. May need to check that VCF and bedgraph are sorted\n";
				$missing_chr{$tok[0]} = 1;
			}
		}
		#update position
		while (!eof($bedfh) && $bedtok[2] < $tok[1]) { bedsplit(\@bedtok, $bedfh); }
	}
	
	## ancestral allele
	if ($alfile) {
		if ($chr ne $tok[0]) {
		# store ancestral alleles in massive matrix
			@outal = ();
			$chr = $1 if ($l1 =~ /^(\S+)/);

			while ($chr ne $tok[0]) {
				my $l1 = <$alfh>;
				die("$tok[0] not found in allele file $alfile\n") if (!$l1);
				$chr = $1 if ($l1 =~ /^(\S+)/);
			}
			$l2 = $l1;
			my $c;

			do {
				push @outal, [(split(/\s+/, $l2))[@col]];
				if (defined($l2 = <$alfh>)) {
					$c = $1 if ($l2 =~ /^(\S+)/);
				}
			} while (defined($l2) && $c eq $chr);

			$l1 = $l2;
			$nsites = scalar(@outal);

			# debug
			#print "$nsites\n";
			#foreach (@outal) {
			#	print "@$_\n";
			#}
			#exit;
		}

		# collect ingroup alleles
		die("More than one reference allele at $tok[0] $tok[1]: $tok[3]\n") if ($tok[3] =~ /,/);
		my $reflen = length($tok[3]);
		my @ingroup = ($tok[3], split(/,/, $tok[4]));

		# collect alleles for each outgroup over the length of the ingroup ref allele
		my @outgroup;
		my $start = $tok[1]-1;
		my $end = $start + $reflen;
		# loop over each outgroup assembly
		for (my $assem=0; $assem<$nalleles; $assem++) {

			if (length($outal[$start]->[$assem]) > $reflen) {
				# insert wrt to ingroup reference
				push @outgroup, $outal[$start]->[$assem];
			} else {
				# non-insertion wrt to ingroup reference
				$outgroup[$assem] = "";
				for (my $s=$start; $s<$end; $s++) {
					$outgroup[$assem] .= $outal[$s]->[$assem] if ($outal[$s]->[$assem] ne '*');
				}

				if (!$outgroup[$assem]) {
					# check if all aligned sites were deleted
					$outgroup[$assem] = '*';
				} else {
					# check for deletions longer that ingroup deletions
					if ($outal[$end-1]->[$assem] eq '*' && $end <= $#outal && $outal[$end]->[$assem] eq '*') {
						$outgroup[$assem] .= '*'; # denotes that next site beyond last reference allele site is deleted
					}
				}
			}
		}

		# find ancestral allele
		my $anc_allele;

		if ($aamethod == 1) {
			# use parsimony to find ancestral allele

			# Find the outgroup matching configuration for each ingroup allele in the lookup table
			my @outsp;
			my $tree;

			foreach my $inallele (@ingroup) {

				my $tabref = \@parstable;
				$tree = "";

				foreach my $outallele (@outgroup) {
					if ($outallele =~ /N/ ) {
						# outgroup allele is missing
						$tabref = ${$tabref}[2];
						$tree .= 2;
					} else {
						# check if ingroup allele matches outgroup
						my $match = $inallele eq $outallele ? 1 : 0;
						$tabref = ${$tabref}[$match];
						$tree .= $match;
					}

				}

				push @outsp, $tabref;
				$tree_counts{$tree}++;
				$tree .= ',';
			}

			# resolve exceptional cases for most parsimonious solution
			if ($nalleles == 3) {
				# 3 outgroup case
				if ($tree =~ /100/ && $tree !~ /011/) {
					push @outsp, 0;
				}
			}

			# pick ancestral allele from the most related outgroup species
			my $spid = (sort {$a <=> $b} @outsp)[0];
			$anc_allele = $spid != 9 ? $outgroup[$spid] : 'N';

			# debug
			#print "$tok[0]\t$tok[1]\t@ingroup;\t@outgroup;\t@outsp;\t$anc_allele\n"; next;			

		} elsif ($aamethod == 2) {
			# use major allele among outgroups as ancestral allele

			my %allele_counts = ();
			foreach my $inallele (@ingroup) {
				$allele_counts{$inallele} = 0;
			}

			# count occurances of outgroup alleles among the ingroup
			my $aa_seen = 0; # records whether at least one ref or alt allele is present among outgroups
			foreach my $outallele (@outgroup) {
				if (exists $allele_counts{$outallele}) {
					$allele_counts{$outallele}++;
					$aa_seen = 1;				
				}
			}

			if ($aa_seen) {
				my $maxcount = 0;
				my $nmax = 0;
					foreach my $allele (@ingroup) {
					if ($allele_counts{$allele} >= $maxcount) {
						$nmax = $allele_counts{$allele} == $maxcount ? $nmax + 1:1;
						$maxcount = $allele_counts{$allele};
						$anc_allele = $allele;
					}
				}
				$anc_allele = "N" if ($nmax > 1 || $maxcount < 1);
		
			} else {
				$anc_allele="N";
			}
		}

		# add ancestral allele to INFO
		if (!($tok[7] =~ s/AA=[^;|\s]+/AA=$anc_allele/)) {
			$tok[7] .= ";AA=$anc_allele";
		}
		$aa_count++ if ($anc_allele ne 'N');
	}

	## genotype respresentation info
	if ($genorep) {
		# reset counts
		foreach my $id (@grpid) {
			$grpcounts{$id}{gq} = 0;
			$grpcounts{$id}{dp} = 0;
		}
		
		# get DP and GQ indices
		my ($dpidx, $gqidx, $iter) = (-9, -9, 0);
		foreach my $fmtid (split(':', $tok[8])) {
			if ($fmtid eq 'DP') {
				$dpidx = $iter;
			} elsif ($fmtid eq 'GQ') {
				$gqidx = $iter;
			}
			last if ($dpidx > -9 && $gqidx > -9);
			$iter++;
		}

		# count genotypes
		for (my $i = 9; $i <= $#tok; $i++) {
			if ($grpmat[$i][0]) {
				my @geno = split(':', $tok[$i]);
				$grpcounts{$grpmat[$i][1]}{dp}++ if ($dpidx > -9 && $geno[0] !~ /\./ && $geno[$dpidx] >= $grpmat[$i][3]);
				$grpcounts{$grpmat[$i][1]}{gq}++ if ($gqidx > -9 && $geno[0] !~ /\./ && $geno[$gqidx] >= $grpmat[$i][2]);
			}
		}

		# record genotype representation info
		my $repdp_str = "REPDP=";
		my $repgq_str = "REPGQ=";
		$iter = 0;
		foreach my $id (@grpid) {
			$repdp_str.="$grpcounts{$id}{dp}" if ($dpidx > -9);
			$repgq_str.="$grpcounts{$id}{gq}" if ($gqidx > -9);
			if ($iter < $#grpid) {
				$repdp_str .= ",";
				$repgq_str .= ",";
			}
			$iter++;
		}
		$repdp_str = "" if ($dpidx == -9);
		$repgq_str = "" if ($gqidx == -9);
		
		addInfo(\$tok[7], "REPDP=", \$repdp_str);
		addInfo(\$tok[7], "REPGQ=", \$repgq_str);
		$tok[7] =~ s/;{2,}/;/g;
		$tok[7] =~ s/;$//;
	}

	## remove filter annotations
	if ($rmfilter) {
		if ($tok[6] ne '.') {
			my $filter_str = "";
			foreach my $fid (split(/;/, $tok[6])) {
				$filter_str .= "$fid;" unless exists $delfilter{$fid};
			}
			$filter_str =~ s/;*$//;
			$tok[6] = (!$filter_str) ? '.' : $filter_str; 
		}
	}

	## apply filters
	if ($dpbounds || $hetbound || $bed) {
		my @filter_annotations;
	
		# coverage filters
		if ($dpbounds) {
			my $dp = $1 if ($tok[7] =~ /DP=(\d+)/);
			if (defined $dp) {
				if ($dp < $dpcutoff[0]) {
					push @filter_annotations, "LowDP";
				} elsif ($dp > $dpcutoff[1]) {
					push @filter_annotations, "HighDP";
				}
			} else {
				print STDERR "No DP info for $tok[0] $tok[1], INFO=<$tok[7]>\n";
			}
		}

		# excess heterozygosity filter
		# Flags a site based on the lowest p-value for multiallelic sites
		if ($hetbound && $tok[7] =~ /ExcHet=([^;]+)/) {
			foreach my $p (split(/,/, $1)) {
				push @filter_annotations, "HighHet" if ($p < $hetbound);
			}
		}

		# bed filters
		if ($bed && $vcfchr eq $bedtok[0] && $tok[1] > $bedtok[1] && $tok[1] <= $bedtok[2]) {
			map {push @filter_annotations, $_ if exists $bedfilter{$_}} split(';',$bedtok[3]); 
		}

		# update FILTER field with new annotations
		if (@filter_annotations) {
			if ($tok[6] eq "." || $tok[6] eq "PASS") {
				$tok[6] = join(';',@filter_annotations);
			} else {
				if ($overwrite) {
					$tok[6] .= join(';',@filter_annotations);
				} else {
					push @filter_annotations, split(';', $tok[6]);
					$tok[6] = join(';', uniq @filter_annotations);
				}
			}
		} else {
			$tok[6] = "PASS" if ($tok[6] eq '.');
		}
	}
	
	if ($tok[6] ne '.' && uc($tok[6]) ne 'PASS') {
		$fail_counts++;
		map {$filtercounts{$_}++} split(';', $tok[6]);
	}
	if (exists $filterconfigs{$tok[6]}) {
		$filterconfigs{$tok[6]}++;
	} else {
		$filterconfigs{$tok[6]} = 1;
	}

	print STDOUT "@tok\n";
}

close $vcffh;
close $alfh if ($alfile);
close $bedfh if ($bed);

# print counts
print STDERR "\n===== REPORT =====\n";
print STDERR "Last contig processed: $tok[0]\n";
print STDERR "Number sites in VCF: $site_count\n";

if ($alfile) {
	print STDERR "\n===== ANCESTRAL ALLELE REPORT =====\n";
	print STDERR "Number of sites with ancestral alleles: $aa_count\n";
	if ($aamethod == 1) {
		print STDERR "Ingroup allele matching configurations\n";
		treeCounter(\%tree_counts, $nalleles, 1);
	}
}

print STDERR "\n===== FILTER SUMMARY =====\n";
print STDERR "Number of FAIL sites: $fail_counts\n";
print STDERR "\n===== FILTER CONFIGURATION COUNTS =====\n";
foreach (keys %filterconfigs) {
	print STDERR "$_: $filterconfigs{$_}\n" unless ($_ eq 'PASS' || $_ eq '.');
}

print STDERR "\n===== FILTER FLAG COUNTS =====\n";
foreach (@{$header{FILTER}}) {
	if ($_ =~ /ID=([^,]+)/) {
		print STDERR "$1: $filtercounts{$1}\n" unless ($1 eq 'PASS' || $1 eq '.');
	}
}

exit;

sub bedsplit {
	my ($tok, $fh) = @_;
	@$tok = split(/\s+/, readline($$fh));
	$$tok[0] = uc($$tok[0]);
}


sub makeLookup {
### returns species representing the ancestral allele ###

### species in descending relatedness to Astatotilapia calliptera
# Simochromis diagramma, Tropheini tribe (0)
# Neolamprologus brichardi, Lamprologini tribe (1)
# Oreochromis niloticus (2)
# not able to assign (9)

### codes corresponding to elements of of the lookup table
# [S. diagramma]->[N. brichardi]->[O. niloticus]
# 0: mismatch to A. calliptera allele
# 1: match A. calliptera allele
# 2: missing data

        my @anc;

        $anc[0]->[0]->[0] = 9;
        $anc[0]->[0]->[1] = 9;
        $anc[0]->[0]->[2] = 9;
        $anc[0]->[1]->[0] = 9;
        $anc[0]->[1]->[1] = 1;
        $anc[0]->[1]->[2] = 1;
        $anc[0]->[2]->[0] = 9;
        $anc[0]->[2]->[1] = 2;
        $anc[0]->[2]->[2] = 9;
        $anc[1]->[0]->[0] = 9;
        $anc[1]->[0]->[1] = 0;
        $anc[1]->[0]->[2] = 0;
        $anc[1]->[1]->[0] = 0;
        $anc[1]->[1]->[1] = 0;
        $anc[1]->[1]->[2] = 0;
        $anc[1]->[2]->[0] = 0;
        $anc[1]->[2]->[1] = 0;
        $anc[1]->[2]->[2] = 0;
        $anc[2]->[0]->[0] = 9;
        $anc[2]->[0]->[1] = 2;
        $anc[2]->[0]->[2] = 9;
        $anc[2]->[1]->[0] = 1;
        $anc[2]->[1]->[1] = 1;
        $anc[2]->[1]->[2] = 1;
        $anc[2]->[2]->[0] = 9;
        $anc[2]->[2]->[1] = 2;
        $anc[2]->[2]->[2] = 9;

        if (@_) {
                my $config_counts = $_[0];
                for (my $i=0; $i<3; $i++) {
                        for (my $j=0; $j<3; $j++) {
                                for (my $k=0; $k<3; $k++) {
                                        $$config_counts{"$i$j$k"} = 0;
                                }
                        }
                }
        }

	return @anc;
}

sub treeCounter {
# fun = 0: initalize count hash
# fun = 1: print count hash

	my ($counts, $nout, $fun) = @_;

	# generate array of configurations
	my @config = ("");

	for (my $i=0; $i<$nout; $i++) {
		my @arr;
		foreach my $str (@config) {
			foreach my $state ((0,1,2)) {
				push @arr, "$str$state";
			}
		}
		@config = @arr;
	}

	# process configurations
	foreach (@config) {
		if ($fun == 0) {
			$$counts{$_} = 0;
		} elsif ($fun == 1) {
			print STDERR "$_: $$counts{$_}\n";
		}
	}
}

sub addInfo {
	my ($info, $oldstr, $newstr) = @_;
	if ($$info =~ /$oldstr/) {
		$$info =~ s/${oldstr}[^;]+/$$newstr/;
	} else {
		$$info .= ";$$newstr";
	}
}
