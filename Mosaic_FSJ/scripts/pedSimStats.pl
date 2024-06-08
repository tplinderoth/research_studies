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

## parse run arguments

die(qq/
pedSimStats.pl

Required inputs:

--ped_file <FILE>      Observed pedigree file.
--indmeta <FILE>       Individual metadata file.
--nestmeta <FILE>      Nest metadata file.
--ancfile <FILE>       File of ancestral individual IDs: column 1 = ID, column 2 = 'resident' or 'translocated'.
--out <STRING>         Output file name prefix.
--simexec <FILE>       simPed.R binary file.
--statsexec <FILE>     relateStats binary file.

Optional inputs:

--rp <FLOAT>           Rate per generation that breeding pairs are replaced by new pairs. [$rp]
--maxoffspring <INT>   Max number of offspring for a breeding pair per generation. [$maxoffspring]
--mature <INT>         Minimum age for an individual to reproduce. [$mature]
--keep_unbanded        Specify to keep unbanded individual in analysis, which are excluded by default.
--p_male <FLOAT>       Probability that individual without known sex is male. [$p_male]
--year <INT>           Year to calculate stats for, default is all years present in ped_file.
--keep_sim             Specify to keep intermediate pedigree simulations, which are deleted by default.
--focal_ind <FILE>     File of focal ancestor IDs \(one per row\) to track contributions for.
--seed_file <FILE>     File of random seeds for simulations \(one per row\). A simulation will be performed for each seed.
--nsims <INT>          Number of simulations to perform. Does not apply if seeds file is provided [$nsims].
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
'nsims=i' => \$nsims);

# check for necessary inputs
die("Must supply --ped_file\n") if (!$ped_file);
die("Must supply --indmeta file\n") if (!$indmeta);
die("Must supply --nestmeta file\n") if (!$nestmeta);
die("Must supply simPed.R binary with --simexec\n") if (!$simexec);
die("Must supply relateStats binary with --statsexec\n") if (!$statsexec);
die("Must specify output file prefix with --out\n") if (!$outprefix);
die("Must supply ancestor IDs with --ancfile\n") if (!$ancfile);

# check validity of necessary inputs
die("ERROR. Unable to find --ped file $ped_file\n") if (!-f $ped_file);
die("ERROR. Unable to find --indmeta file $indmeta\n") if (!-f $indmeta);
die("ERROR. Unable to find --nestmeta file $nestmeta\n") if (!-f $nestmeta);
die("ERROR. Unable to find --ancfile file $ancfile\n") if (!-f $ancfile);
die("ERROR. Unable to find simPed.R binary\n") if (!-e $simexec);
die("ERROR. simPed.R binary is not executable") if (!-x $simexec);
die("ERROR. Unable to find relateStats binary\n") if (!-e $statsexec);
die("ERROR. relateStats is not executable\n") if (!-x $statsexec);

# check optional args
die("ERROR. --rp $rp out of range, must be in [0,1]\n") if ($rp < 0 || $rp > 1);
die("ERROR. Invalid --maxoffspring $maxoffspring, must be a positive integer\n") if ($maxoffspring < 0);
die("ERROR. Invalid --mature $mature, must be a positive integer\n") if ($mature < 0);
die("ERROR. --p_male $p_male out of range, must be in [0,1]\n") if ($p_male < 0 || $p_male > 1);
die("ERROR. Invalid --year $year, must be a positive integer\n") if ($year && $year < 0);
die("ERROR. Invalid --nsims $nsims, must be a positive integer\n") if ($nsims < 0);
die("ERROR. Unable to find --focal_ind file $focal_ind\n") if ($focal_ind && !-f $focal_ind);
die("ERROR. --focal_ind requires --ancfile\n") if ($focal_ind && !$ancfile);

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

## initialize array to hold stats

my @stats = (0) x 3; # [max_contribution, n_resident_contributors, n_trans_contributors, focal individual contributions (in input order)]

## add slots in stats array for focal individual values
my @fid;
my %focal_hash;
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

## set up output streams
my @outfh;
foreach my $i (0 .. $#years) {
	my $outname = "${outprefix}_${years[$i]}.tsv";
	open($outfh[$i], '>>', "$outname") or die("Unable to open output filehandle $outname\n");
	print { $outfh[$i] } "SEED\tMAX_CONTRIBUTION\tN_RESIDENT_CONTRIBUTORS\tN_TRANSLOCATED_CONTRIBUTORS\t@fid\n"; # header
}

## simulate pedigrees and calculate stats
my $simopts = "--rp $rp --maxoffspring $maxoffspring --mature $mature --p_male $p_male --rmatrix 1";
$simopts .= " --keep_unbanded" if ($keep_unbanded);
$simopts .= " --anc $ancfile" if ($ancfile);
my %anc_hash;
my %sim_focal_hash;

foreach my $seed (@seeds) {
	# simulate pedigree
	my $simout = "${outprefix}_${seed}";
	my $simcmd = "$simexec $simopts --seed $seed --out $simout $ped_file $indmeta $nestmeta";
	#print STDERR "$simcmd\n";
	my $rv = system($simcmd);
	die("ERROR. Failure executing pedigree simulation: $simcmd\n") if ($rv);

	# calculate contributions to population over time
	my $out_file = 0;
	foreach my $yr (@years) {
		my $sim_anc = "${simout}.id.anc";
		system("cut -f1 ${simout}.id | tail -n+2 > $sim_anc");
		my $statout = "${simout}_${yr}";
		my $statcmd = "$statsexec --pedstat 1 --ped ${simout}.ped --rmat ${simout}.mat --anc $sim_anc --time2 $yr --out $statout";
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

		# update iterators		
		$out_file++;
	}

	# delete simulated pedigree files
	unless ($keep_sims) {
		unlink("${simout}.id");
		unlink("${simout}.id.anc");
		unlink("${simout}.mat");
		unlink("${simout}.ped");
	}
	
}

foreach my $i (0 .. $#years) { 
	close $outfh[$i]; 
}

exit 0;

## define subroutines

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
