#!/usr/bin/perl

# pedSimStats.pl

use warnings;
use strict;
use Getopt::Long;
use List::MoreUtils qw(uniq);

## specify parameters

my $seed_file = undef;
my $ped_file = undef;
my $indmeta = undef;
my $nestmeta = undef;
my $rp = 0;
my $ancfile = undef;
my $maxoffspring = 5;
my $mature = 2;
my $keep_unbanded = undef;
my $p_male = 0.5;
my $year = undef;
my $keep_sims = 0;
my $outprefix = undef;
my $focal_ind = undef;
my $simexec = undef;
my $statsexec = undef;
my $nsims = 1;
my $cohort = 0;
my $doContribution = undef;
my $doInbreed = undef;
my $gp_prune = 1;

## parse run arguments

die(qq/
pedSimStats.pl

Universal required inputs:

--ped_file <FILE>      Observed pedigree file.
--indmeta <FILE>       Individual metadata file.
--nestmeta <FILE>      Nest metadata file.
--out <STRING>         Output file name prefix.
--simexec <FILE>       simPed.R binary file.


Analysis types:

--doContribution      Expected genetic contributions.
Requires:
--statsexec <FILE>     relateStats binary file.
--ancfile <FILE>       File of ancestral individual IDs: column 1 = ID, column 2 = 'resident' or 'translocated'.
Optional:
--focal_ind <FILE>     File of focal ancestor IDs \(one per row\) to track contributions for.
--cohort               Restrict ancestor contributions to local recruits born in the respective year.

--doInbreed            Average pedigree inbreeding among local recurit cohorts and population.
Optional:
--gp_prune <0|1>       Calculate average F using 1 = only individuals with all known grandparents, or 0 = all individuals. [$gp_prune] 
Note: Input ped file can have 'F' column, which will set founder inbreeding.

Other optional inputs:

--rp <FLOAT>           Rate per generation that breeding pairs are replaced by new pairs. [$rp]
--maxoffspring <INT>   Max number of offspring for a breeding pair per generation. [$maxoffspring]
--mature <INT>         Minimum age for an individual to reproduce. [$mature]
--keep_unbanded        Specify to keep unbanded individual in analysis, which are excluded by default.
--p_male <FLOAT>       Probability that individual without known sex is male. [$p_male]
--year <INT>           Year to calculate stats for, default is all years present in ped_file.
--keep_sim             Specify to keep intermediate pedigree simulations, which are deleted by default.
--seed_file <FILE>     File of random seeds for simulations \(one per row\). A simulation will be performed for each seed.
--nsims <INT>          Number of simulations to perform. Does not apply if seeds file is provided [$nsims].


If no analysis type is provided this script can be used to just generate pedigrees with --keep_sim
\n/) if (!@ARGV);

GetOptions(
'ped_file=s' => \$ped_file,
'indmeta=s' => \$indmeta,
'nestmeta=s' => \$nestmeta,
'out=s' => \$outprefix,
'simexec=s' => \$simexec,
'statsexec=s' => \$statsexec,
'rp=f' => \$rp,
'ancfile=s' => \$ancfile,
'maxoffspring=i' => \$maxoffspring,
'mature=i' => \$mature,
'keep_unbanded' => \$keep_unbanded,
'p_male=f' => \$p_male,
'year=i' => \$year,
'keep_sims' => \$keep_sims,
'focal_ind=s' => \$focal_ind,
'seed_file=s' => \$seed_file,
'nsims=i' => \$nsims,
'cohort' => \$cohort,
'gp_prune=i' => \$gp_prune,
'doContribution' => \$doContribution,
'doInbreed' => \$doInbreed);

# check for universal necessary inputs
die("Must supply --ped_file\n") if (!$ped_file);
die("Must supply --indmeta file\n") if (!$indmeta);
die("Must supply --nestmeta file\n") if (!$nestmeta);
die("Must supply simPed.R binary with --simexec\n") if (!$simexec);
die("Must specify output file prefix with --out\n") if (!$outprefix);

# check validity of necessary inputs
die("ERROR. Unable to find --ped file $ped_file\n") if (!-f $ped_file);
die("ERROR. Unable to find --indmeta file $indmeta\n") if (!-f $indmeta);
die("ERROR. Unable to find --nestmeta file $nestmeta\n") if (!-f $nestmeta);
die("ERROR. Unable to find simPed.R binary\n") if (!-e $simexec);
die("ERROR. simPed.R binary is not executable") if (!-x $simexec);

# check contribution calculation inputs
if ($doContribution) {
	# necessary
	die("--doContribution requires relateStats binary with --statsexec\n") if (!$statsexec);
	die("ERROR. Unable to find relateStats binary\n") if (!-e $statsexec);
	die("ERROR. relateStats is not executable\n") if (!-x $statsexec);
	die("--doContribution requires ancestor IDs with --ancfile\n") if (!$ancfile);
	die("ERROR. Unable to find --ancfile file $ancfile\n") if (!-f $ancfile);

	# optional
	die("ERROR. Unable to find --focal_ind file $focal_ind\n") if ($focal_ind && !-f $focal_ind);
	die("ERROR. --focal_ind requires --ancfile\n") if ($focal_ind && !$ancfile);
}

# check optional args
die("ERROR. --rp $rp out of range, must be in [0,1]\n") if ($rp < 0 || $rp > 1);
die("ERROR. Invalid --maxoffspring $maxoffspring, must be a positive integer\n") if ($maxoffspring < 0);
die("ERROR. Invalid --mature $mature, must be a positive integer\n") if ($mature < 0);
die("ERROR. --p_male $p_male out of range, must be in [0,1]\n") if ($p_male < 0 || $p_male > 1);
die("ERROR. Invalid --year $year, must be a positive integer\n") if ($year && $year < 0);
die("ERROR. Invalid --nsims $nsims, must be a positive integer\n") if ($nsims < 0);
die("ERROR. --gp_prune must be either 0 or 1") unless ($gp_prune == 0 || $gp_prune == 1);

## set style for printing arrays
$" = "\t";

## prepare seeds
my @seeds;

my $seedfh = undef;
if ($seed_file) {
	open($seedfh, '<', $seed_file) or die("Unable to open seed file $seed_file: $!\n");
}

if ($seedfh) {
	chomp(@seeds = <$seedfh>);
	close $seedfh;
} else {
	for (my $i=0; $i < $nsims; $i++) {
		@seeds = map {1 + int(rand($nsims + 1e6))} 1..$nsims; 
	}
}

close $seedfh;


## determine years to analyze
my @years;
if ($year) {
	push @years, $year;
} else {
	open(my $pedfh, '<', $ped_file) or die("Unable to open ped file $ped_file: $!\n");
	chomp(my $head = <$pedfh>);
	my @headarr = split(/\s+/, $head);
	my $cohort_idx = 0;
	foreach (@headarr) {
		last if (uc($_) eq "COHORT");
		$cohort_idx++;
	}
	die("COHORT field not found in ped file\n") if ($cohort_idx > $#headarr);
	while (<$pedfh>) {
		chomp;
		my @tok = split(/\s+/, $_);
		push @years, $tok[$cohort_idx];
	}

	@years = sort {$a <=> $b} uniq(@years);
}

## initialize array to hold average inbreeding values and inbreeding output filehandles
my @fcohort = ('NaN') x scalar(@years);
my @fpop = @fcohort;
my @ffh;
my (%f_cohort, %f_pop);
if ($doInbreed) {
	my $outf_name = "${outprefix}.fcohort";
	open($ffh[0], '>>', $outf_name) or die("Unable to open output filehandle $outf_name\n");
	$outf_name = "${outprefix}.fpop";
	open($ffh[1], '>>', $outf_name) or die("Unable to open output filehandle $outf_name\n");
	print { $ffh[0] } "SEED\tF", join("\tF", @years), "\n";
	print { $ffh[1] } "SEED\tF", join("\tF", @years), "\n";
	foreach my $yr (@years) {
		@{$f_cohort{$yr}} = ();
		@{$f_pop{$yr}} = ();
	}
}

## initialize array to hold contribution stats and contribution output filehandles
my @stats = (0) x 3; # [max_contribution, n_resident_contributors, n_trans_contributors, focal individual contributions (in input order)] - these are contribution stats
my @outfh;
my %focal_hash;

if ($doContribution) {
	my @fid;
	# add slots in stats array for focal individual values	
	if ($focal_ind) {
		open(my $focalfh, '<', $focal_ind) or die("Unable to open file of focal individuals $focal_ind: $!\n");
		chomp(@fid = <$focalfh>);
		close $focalfh;
		my $i = $#stats + 1;
		foreach (@fid) {
			push @stats, 0;
			$focal_hash{$_} = $i;
			$i++;
		}
	}

	## set up output streams
	foreach my $i (0 .. $#years) {
        	my $outname = "${outprefix}_${years[$i]}_contribution.tsv";
        	open($outfh[$i], '>>', "$outname") or die("Unable to open output filehandle $outname\n");
        	print { $outfh[$i] } "SEED\tMAX_CONTRIBUTION\tN_RESIDENT_CONTRIBUTORS\tN_TRANSLOCATED_CONTRIBUTORS\t@fid\n"; # header
	}
}

## simulate pedigrees and calculate stats
my $simopts = "--rp $rp --maxoffspring $maxoffspring --mature $mature --p_male $p_male";
$simopts .= " --rmatrix 1" if ($doContribution);
$simopts .= " --calc_f" if ($doInbreed);
$simopts .= " --keep_unbanded" if ($keep_unbanded);
$simopts .= " --anc $ancfile" if ($ancfile);
my $simped_fh;
my %anc_hash;
my %sim_focal_hash;
my %ped;
my $fidx = undef;

foreach my $seed (@seeds) {
	# simulate pedigree
	my $simout = "${outprefix}_${seed}";
	my $simcmd = "$simexec $simopts --seed $seed --out $simout $ped_file $indmeta $nestmeta";
	#print STDERR "$simcmd\n";
	my $rv = system($simcmd);
	die("ERROR. Failure executing pedigree simulation: $simcmd\n") if ($rv);

	# calculate contributions to population over time
	if ($doContribution) {
		my $out_file = 0;
		foreach my $yr (@years) {
			# Generate list of local recruits born in present year if analyzing cohort
			my $cohort_f;
			my $n_lr = 0;
			if ($cohort) {
				open($simped_fh, '<', "${simout}.ped") or die("ERROR. Unable to open pedigree simulation file ${simout}.ped");
				$cohort_f = "${simout}.cohort";
				open(my $cohort_fh, '>', $cohort_f) or die("ERROR. Unable to open file to store cohort IDs: $cohort_f\n");
				<$simped_fh>; # skip header
				while (my $pedline = <$simped_fh>) {
					chomp($pedline);
					my @tok = split(/\t/, $pedline);
					if ($tok[4] == $yr && $tok[6] eq "CORE_LR") {
						print $cohort_fh "$tok[0]\n";
						$n_lr++;
					}
				}
				close $simped_fh;
				close $cohort_fh;
			}

			my $sim_anc = "${simout}.id.anc";
			system("cut -f1 ${simout}.id | tail -n+2 > $sim_anc");
			my $statout = "${simout}_${yr}";
			my $statcmd = "$statsexec --pedstat 1 --ped ${simout}.ped --rmat ${simout}.mat --anc $sim_anc --out $statout";
			$statcmd .= $cohort ? " --cohort $cohort_f" : " --time2 $yr";

			if ($cohort && $n_lr < 1) {
				# no local recruits to calculate contributions to
				foreach my $st (@stats) {$st = "NA";}
				print { $outfh[$out_file] } "$seed\t@stats\n";
				$out_file++;
				next;
			}

			#print STDERR "$statcmd\n";
			$rv = system($statcmd);
			die("ERROR. Failure calculating genetic contributions: $statcmd\n") if ($rv);

			# add ancestor origin to hash for quick checking and map stats index to individual's sim ID
			%anc_hash = ();
			%sim_focal_hash = ();
			open(my $idfh, '<', "${simout}.id") or die("Unable to open ID map file ${simout}.id\n");
			<$idfh>; # skip header
			while (my $idline = <$idfh>) {
				chomp $idline;
				my @tok = split(/\t/, $idline);
				my $origin = uc($tok[2]);
				if ($origin eq "CORE_RESIDENT") {
					$anc_hash{$tok[0]} = 1;
				} elsif ($origin eq "SITE_1" || $origin eq "SITE_12" || $origin eq "SITE_13" || $origin eq "SITE_18" || $origin eq "TEXACO") {
					$anc_hash{$tok[0]} = 2;
				} else {
					$anc_hash{$tok[0]} = 3;
				}
				$sim_focal_hash{$tok[0]} = $focal_hash{$tok[1]} if ($focal_ind && exists $focal_hash{$tok[1]});
			}
			close $idfh;

			# calculate contribution stats
			open(my $cfh, '<', "${statout}.pedstat1") or die("Unable to open pedstat file ${statout}.pedstat1\n");
			contributionStats($cfh, \@stats, \%anc_hash, \%sim_focal_hash);

			# print stats to output
			print { $outfh[$out_file] } "$seed\t@stats\n";

			# delete relateStats files
			unless ($keep_sims) {
				unlink("${statout}.pedstat1");
			}
			unlink($cohort_f) if ($cohort && -f $cohort_f);

			# update iterators		
			$out_file++;
		}
	}

	if ($doInbreed) {
		open($simped_fh, '<', "${simout}.ped") or die("ERROR. Unable to open pedigree simulation file ${simout}.ped");
		collectF($simped_fh, \%f_cohort, \%f_pop, \$fidx, \%ped, $gp_prune);
		# calculate cohort means
		meanF(\@fcohort, \%f_cohort, \@years);
		# calculate population means
		meanF(\@fpop, \%f_pop, \@years);
		# print output
		print { $ffh[0] } "$seed\t@fcohort\n";
		print { $ffh[1] } "$seed\t@fpop\n";
	}

	# delete simulated pedigree files
	unless ($keep_sims) {
		unlink("${simout}.id") if (-f "${simout}.id");
		unlink("${simout}.id.anc") if (-f "${simout}.id.anc");
		unlink("${simout}.mat") if (-f "${simout}.mat");
		unlink("${simout}.ped");
	}
	
}

# close output filehandles

if (@outfh) {
	foreach my $i (0 .. $#years) { 
		close $outfh[$i]; 
	}
}

exit 0;

## define subroutines

sub meanF {
	my ($fmean, $fvals, $years) = @_;
	
	foreach my $val (@{$fmean}) {$val = 'NaN';}

	my $i = 0;
	foreach my $yr (@{$years}) {
		if (exists $$fvals{$yr}) {
			my $n = scalar(@{$$fvals{$yr}});
			if ($n > 0) {
				my $sum = 0.0;
				foreach (@{$$fvals{$yr}}) {$sum += $_;}
				$$fmean[$i] = $sum/$n;
			}
		}
		$i++;
	}
}

sub collectF {
	my ($fh, $f_cohort, $f_pop, $fidx, $ped, $reqgp) = @_;

	foreach (keys %{$f_cohort}) {
		@{$$f_cohort{$_}} = ();
	}
	foreach (keys %{$f_pop}) {
		@{$$f_pop{$_}} = ();
	}

	# collect index of pedigree column if not already known
	my $found = 0;
	if (!$$fidx) {
		$$fidx = 0;
		chomp(my $h = <$fh>);
		foreach my $field (split(/\t/,$h)) {
			if  (uc($field) eq 'PEDIGREE_F') {
				$found = 1;
				last;
			}
			$$fidx++;
		}
	} else {
		$found = 1;
		<$fh>; # skip header
	}
	die("ERROR in collectF(): No column 'PEDIGREE_F' column in input ped\n") if (!$found);

	# store pedigree
	%$ped = ();
	while (<$fh>) {
		chomp;
		my @tok = split(/\t/, $_);
		next if ($tok[0] eq '*'); # there should not be missing main IDs, but just in case
		$$ped{$tok[0]}{fid} = $tok[1];
		$$ped{$tok[0]}{mid} = $tok[2];
		$$ped{$tok[0]}{cohort} = $tok[4];
		$$ped{$tok[0]}{cohort_last} = $tok[5];
		$$ped{$tok[0]}{pop} = $tok[6];
		$$ped{$tok[0]}{f} = $tok[$$fidx];
	}

	# store F coefficient for Core Region local recruits for which both grandparents are known
	foreach my $id (keys %$ped) {
		if ($$ped{$id}{pop} eq 'CORE_LR') {
			if ($reqgp) {
				my $father = $$ped{$id}{fid};
				my $mother = $$ped{$id}{mid};
				if (exists $$ped{$father} && exists $$ped{$mother}) {
					# check that both granparents exist
					next unless (exists $$ped{$$ped{$father}{fid}} && exists $$ped{$$ped{$father}{mid}} && exists $$ped{$$ped{$mother}{fid}} && exists $$ped{$$ped{$mother}{mid}});
				} else {
					next;
				}
			}
			push @{$f_cohort{$$ped{$id}{cohort}}}, $$ped{$id}{f};
			foreach my $yr ($$ped{$id}{cohort} .. $$ped{$id}{cohort_last}) {
				push @{$f_pop{$yr}}, $$ped{$id}{f} if exists $$f_pop{$yr};
			}
		}
	}
}

sub contributionStats {
	my ($fh, $stats, $anc, $focal) = @_;

	foreach my $val (@{$stats}) {$val = 0;}

	$$stats[0] = -1; # max contribution

	<$fh>; # skip header
	while (<$fh>) {
		chomp;
		my @tok = split(/\t/,$_);
		# check for max
		$$stats[0] = $tok[2] if ($tok[2] =~ /\d/ && $tok[2] > $$stats[0]);
		# iterate conributors
		$$stats[1]++ if ($tok[2] > 0 && exists $$anc{$tok[0]} && $$anc{$tok[0]} == 1);
		$$stats[2]++ if ($tok[2] > 0 && exists $$anc{$tok[0]} && $$anc{$tok[0]} == 2);
		# log focal individual contribution
		$$stats[$$focal{$tok[0]}] = $tok[2] if (exists $$focal{$tok[0]});
	}
}
