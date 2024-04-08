#!/usr/bin/perl

# ped_problems.pl

use warnings;
use strict;

my $survivorship_wnesting = '/home/tyler/Desktop/tmp/fsj/ped_full.csv'; # from M4_FSJ_Blood_Samples.xlsx SURVIVORSHIP_WNESTING tab
my $ped_source = '/home/tyler/Desktop/tmp/fsj/ped_source.csv'; # from 'JAY DATA WITH PARENTS through 2021.xlsm' JAYID tab
my $pedigree = '/home/tyler/Dropbox/research/jay/mosaic/pedigree/mosaic_ped_contribution_input.txt'; # input for pedigree analysis

open(my $snfh, '<', $survivorship_wnesting) or die("Unable to open file $survivorship_wnesting: $!\n");
open(my $sourcefh, '<', $ped_source) or die("Unable to open file $ped_source: $!\n");
open(my $pedfh, '<', $pedigree) or die("Unable to open file $pedigree: $!\n");

my %data;

# process pedigree source file

<$sourcefh>; # skip header
while (<$sourcefh>) {
	chomp;
	my @tok = split('\t', $_); 
	my $id = $tok[0] ? lc($tok[0]) : "";
	my $mom = (!$tok[3] || $tok[3] =~ /#N\/A/g) ? "" : lc($tok[3]);
	my $dad = (!$tok[4] || $tok[4] =~ /#N\/A/g) ? "" : lc($tok[4]);

	if ($id) {
		if (!exists $data{$id}{source}) {
			@{$data{$id}{source}{sireof}} = ();
			@{$data{$id}{source}{damof}} = ();
		}
		$data{$id}{source}{sex} = $tok[2];
		$data{$id}{source}{mother} = $mom;
		$data{$id}{source}{father} = $dad;
	}

	if ($mom && $mom ne "unknown") {
		if (!exists $data{$mom}{source}) {
			$data{$mom}{source}{sex} = "";
			$data{$mom}{source}{mother} = "";
			$data{$mom}{source}{father} = "";
			@{$data{$mom}{source}{damof}} = ();
			@{$data{$mom}{source}{sireof}} = ();
		}
		push @{$data{$mom}{source}{damof}}, $id if ($id);
	}

	if ($dad && $dad ne "unknown") {
		if (!exists $data{$dad}{source}) {
			$data{$dad}{source}{sex} = "";
			$data{$dad}{source}{mother} = "";
			$data{$dad}{source}{father} = "";
			@{$data{$dad}{source}{damof}} = ();
			@{$data{$dad}{source}{sireof}} = ();
		}
		push @{$data{$dad}{source}{sireof}}, $id if ($id);
	}

}
close $sourcefh;

# process survivorship w/ nesting file

<$snfh>; # skip header
while (<$snfh>) {
	chomp;
	my @tok = split('\t', $_);
	my $id = $tok[0] ? lc($tok[0]) : "";
	my $mom = (!$tok[3]) ? "" : lc($tok[3]);
	my $dad = (!$tok[4]) ? "" : lc($tok[4]);

	if ($id) {
		if (!exists $data{$id}{survive_nest}) {
			@{$data{$id}{survive_nest}{sireof}} = ();
			@{$data{$id}{survive_nest}{damof}} = ();
		}
		$data{$id}{survive_nest}{sex} = $tok[2];
		$data{$id}{survive_nest}{mother} = $mom;
		$data{$id}{survive_nest}{father} = $dad;
	}

	if ($mom && $mom ne "unknown") {
		if (!exists $data{$mom}{survive_nest}) {
			$data{$mom}{survive_nest}{sex} = "";
			$data{$mom}{survive_nest}{mother} = "";
			$data{$mom}{survive_nest}{father} = "";
			@{$data{$mom}{survive_nest}{damof}} = ();
			@{$data{$mom}{survive_nest}{sireof}} = ();
		}
		push @{$data{$mom}{survive_nest}{damof}}, $id if ($id);
	}

	if ($dad && $dad ne "unknown") {
		if (!exists $data{$dad}{survive_nest}) {
			$data{$dad}{survive_nest}{sex} = "";
			$data{$dad}{survive_nest}{mother} = "";
			$data{$dad}{survive_nest}{father} = "";
			@{$data{$dad}{survive_nest}{damof}} = ();
			@{$data{$dad}{survive_nest}{sireof}} = ();
		}
		push @{$data{$dad}{survive_nest}{sireof}}, $id if ($id);
	}

}

close $snfh;

# process pedigree input file

<$pedfh>; # skip header
while (<$pedfh>) {
	chomp;
	my @tok = split('\t', $_);
	my $id = ($tok[0] && $tok[0] ne "0")? lc($tok[0]) : "";
	my $mom = (!$tok[1] || $tok[1] eq "0") ? "" : lc($tok[1]);
	my $dad = (!$tok[2] || $tok[2] eq "0") ? "" : lc($tok[2]);

	if ($id) {
		if (!exists $data{$id}{pedin}) {
			@{$data{$id}{pedin}{sireof}} = ();
			@{$data{$id}{pedin}{damof}} = ();
		}
		$data{$id}{pedin}{sex} = "";
		if ($tok[3] && $tok[3] eq "1") {
			$data{$id}{pedin}{sex} = "male";
		} elsif ($tok[3] && $tok[3] eq "2") {
			$data{$id}{pedin}{sex} = "female";		
		}
		$data{$id}{pedin}{mother} = $mom;
		$data{$id}{pedin}{father} = $dad;
	}

	if ($mom && $mom ne "0") {
		if (!exists $data{$mom}{pedin}) {
			$data{$mom}{pedin}{sex} = "";
			$data{$mom}{pedin}{mother} = "";
			$data{$mom}{pedin}{father} = "";
			@{$data{$mom}{pedin}{damof}} = ();
			@{$data{$mom}{pedin}{sireof}} = ();
		}
		push @{$data{$mom}{pedin}{damof}}, $id if ($id);
	}

	if ($dad && $dad ne "0") {
		if (!exists $data{$dad}{pedin}) {
			$data{$dad}{pedin}{sex} = "";
			$data{$dad}{pedin}{mother} = "";
			$data{$dad}{pedin}{father} = "";
			@{$data{$dad}{pedin}{damof}} = ();
			@{$data{$dad}{pedin}{sireof}} = ();
		}
		push @{$data{$dad}{pedin}{sireof}}, $id if ($id);
	}
}

close $pedfh;

# check input debug
#foreach my $id (keys %data) {
	#if (exists $data{$id}{source}) {
	#	print STDOUT "$id\tsex:$data{$id}{source}{sex}\tmother:$data{$id}{source}{mother}\tfather:$data{$id}{source}{father}\tdamof:@{$data{$id}{source}{damof}}\tsireof:@{$data{$id}{source}{sireof}}\n";
	#} else {
	#	print STDOUT "$id does not exist in SOURCE\n";	
	#}
#	if (exists $data{$id}{survive_nest}) {
#		print STDOUT "$id\tsex:$data{$id}{survive_nest}{sex}\tmother:$data{$id}{survive_nest}{mother}\tfather:$data{$id}{survive_nest}{father}\tdamof:@{$data{$id}{survive_nest}{damof}}\tsireof:@{$data{$id}{survive_nest}{sireof}}\n";
#	} else {
#		print STDOUT "$id does not exist in SURVIVORSHIP_WNESTING\n";
#	}

#	if (exists $data{$id}{pedin}) {
#		print STDOUT "$id\tsex:$data{$id}{pedin}{sex}\tmother:$data{$id}{pedin}{mother}\tfather:$data{$id}{pedin}{father}\tdamof:@{$data{$id}{pedin}{damof}}\tsireof:@{$data{$id}{pedin}{sireof}}\n";
#	} else {
#		print STDOUT "$id does not exist in PEDIGREE INPUT\n";
#	}
#}


### Look for problems and discrepancies

## both dam and sire
my $header = 0;
foreach my $id (keys %data) {
	my $err = 0;
	my $errstr = "";
	if (exists $data{$id}{source}) {
		$err++ if (@{$data{$id}{source}{damof}} && @{$data{$id}{source}{sireof}});
		$errstr .= " [SOURCE: dam_of: @{$data{$id}{source}{damof}} sire_of: @{$data{$id}{source}{sireof}}]";
	}
	if (exists $data{$id}{pedin}) {
		$err++ if (@{$data{$id}{pedin}{damof}} && @{$data{$id}{pedin}{sireof}});
		$errstr .= " [PED: dam_of: @{$data{$id}{pedin}{damof}} sire_of: @{$data{$id}{pedin}{sireof}}]";
	}
	if (exists $data{$id}{survive_nest}) {
		$err++ if (@{$data{$id}{survive_nest}{damof}} && @{$data{$id}{survive_nest}{sireof}});
		$errstr .= " [SURVIVORSHIP_WNESTING: dam_of: @{$data{$id}{survive_nest}{damof}} sire_of: @{$data{$id}{survive_nest}{sireof}}]";
	}

	if ($err) {
		if (!$header) {
			print STDOUT "### BOTH DAM AND SIRE ###\n";
			$header++;
		}
		$errstr =~ s/^ //;
		print STDOUT "$id $errstr\n";
	}
}

## different offspring
$header = 0;
foreach my $id (keys %data) {
	my $err = 0;
	my $errstr = "";
	my $source_offspring = "";
	my $ped_offspring = "";
	my $sn_offspring = "";
	if (exists $data{$id}{source}) {
		$source_offspring = join('_', sort {$a cmp $b} @{$data{$id}{source}{damof}}) . "_"  . join('_', sort {$a cmp $b} @{$data{$id}{source}{sireof}});
	}
	if (exists $data{$id}{pedin}) {
		$ped_offspring = join('_', sort {$a cmp $b} @{$data{$id}{pedin}{damof}}) . "_"  . join('_', sort {$a cmp $b} @{$data{$id}{pedin}{sireof}});
	}
	if (exists $data{$id}{survive_nest}) {
		$sn_offspring = join('_', sort {$a cmp $b} @{$data{$id}{survive_nest}{damof}}) . "_"  . join('_', sort {$a cmp $b} @{$data{$id}{survive_nest}{sireof}});
	}

	if ($source_offspring ne $ped_offspring || $source_offspring ne $sn_offspring || $ped_offspring ne $sn_offspring) {
		$err++;
		if (exists $data{$id}{source}) {
			my @d = sort {$a cmp $b} @{$data{$id}{source}{damof}};
			my @s = sort {$a cmp $b} @{$data{$id}{source}{sireof}};
			$errstr .= " [SOURCE: dam_of: @d sire_of: @s]";
		}
		if (exists $data{$id}{pedin}) {
			my @d = sort {$a cmp $b} @{$data{$id}{pedin}{damof}};
			my @s = sort {$a cmp $b} @{$data{$id}{pedin}{sireof}};
			$errstr .= " [PED: dam_of: @d sire_of: @s]";
		}
		if (exists $data{$id}{survive_nest}) {
			my @d = sort {$a cmp $b} @{$data{$id}{survive_nest}{damof}};
			my @s = sort {$a cmp $b} @{$data{$id}{survive_nest}{sireof}};
			$errstr .= " [SURVIVORSHIP_WNESTING: dam_of: @d sire_of: @s]";
		}
	}

	if ($err) {
		if (!$header) {
			print STDOUT "\n\n### DIFFERENT OFFSPRING ###\n";
			$header++;
		}
		$errstr =~ s/^ //;
		print STDOUT "$id $errstr\n";
	}
}

## different parents
$header = 0;
foreach my $id (keys %data) {
	my $err = 0;
	my $errstr = "";
	my $source_parents = exists $data{$id}{source} ? $data{$id}{source}{mother} . "_" . $data{$id}{source}{father} : "";
	my $ped_parents = exists $data{$id}{pedin} ? $data{$id}{pedin}{mother} . "_" . $data{$id}{pedin}{father} : "";
	my $sn_parents = exists $data{$id}{survive_nest} ? $data{$id}{survive_nest}{mother} . "_" . $data{$id}{survive_nest}{father} : "";;
	
	if ($source_parents ne $ped_parents || $source_parents ne $sn_parents || $ped_parents ne $sn_parents) {
		$err++;
		$errstr .= " [SOURCE: mother: $data{$id}{source}{mother} father: $data{$id}{source}{father}]" if exists $data{$id}{source};
		$errstr .= " [PED: mother: $data{$id}{pedin}{mother} father: $data{$id}{pedin}{father}]" if exists $data{$id}{pedin};
		$errstr .= " [SURVIVORSHIP_WNESTING: mother: $data{$id}{survive_nest}{mother} father: $data{$id}{survive_nest}{father}]" if exists $data{$id}{survive_nest};
	}

	if ($err) {
		if (!$header) {
			print STDOUT "\n\n### DIFFERENT PARENTS ###\n";
			$header++;
		}
		$errstr =~ s/^ //;
		print STDOUT "$id $errstr\n";
	}
}


## different sex
$header = 0;
foreach my $id (keys %data) {
	my $err = 0;
	my $errstr = "";

	my $source_sex = "NA";
	my $ped_sex = "NA";
	my $sn_sex = "NA";

	if (exists $data{$id}{source}) {
		$source_sex = lc($data{$id}{source}{sex});
		$source_sex = "NA" if ($source_sex eq "" || $source_sex eq "unknown");	
	}
	if (exists $data{$id}{pedin}) {
		$ped_sex = lc($data{$id}{pedin}{sex});
		$ped_sex = "NA" if ($ped_sex eq "" || $ped_sex eq "unknown");	
	}
	if (exists $data{$id}{survive_nest}) {
		$sn_sex = lc($data{$id}{survive_nest}{sex});
		$sn_sex = "NA" if ($sn_sex eq "" || $sn_sex eq "unknown");	
	}

	if ((exists $data{$id}{source} && exists $data{$id}{pedin} && $source_sex ne $ped_sex) || (exists $data{$id}{source} && exists $data{$id}{survive_nest} && $source_sex ne $sn_sex) || (exists $data{$id}{pedin} && exists $data{$id}{survive_nest} && $ped_sex ne $sn_sex)) {
		$err++;
		$errstr .= " [SOURCE sex: $source_sex]" if exists $data{$id}{source};
		$errstr .= " [PED sex: $ped_sex]" if exists $data{$id}{pedin};
		$errstr .= " [SURVIVORSHIP_WNESTING sex: $sn_sex]" if exists $data{$id}{survive_nest};
	}

	if ($err) {
		if (!$header) {
			print STDOUT "\n\n### DIFFERENT SEX ###\n";
			$header++;
		}
		$errstr =~ s/^ //;
		print STDOUT "$id $errstr\n";
	}
}

## missing individuals
$header = 0;
foreach my $id (keys %data) {
	if (!exists $data{$id}{source}) {
		print STDOUT "\n\n### MISSING FROM SOURCE (file: 'JAY ID DATA WITH PARENTS through 2021.xlsm', tab: 'JAYID') ###\n" if !$header;
		$header++;
		print STDOUT "$id\n";
	}
}

$header = 0;
foreach my $id (keys %data) {
	if (!exists $data{$id}{pedin}) {
		print STDOUT "\n\n### MISSING FROM PEDIGREE INPUT (file: 'mosaic_ped.csv') ###\n" if !$header;
		$header++;
		print STDOUT "$id\n";
	}
}

$header = 0;
foreach my $id (keys %data) {
	if (!exists $data{$id}{survive_nest}) {
		print STDOUT "\n\n### MISSING FROM SURVIVORSHIP_WNESTING METADATA (file: 'M4_FSJ_Blood_Samples.xlsx', tab: 'SURVIVORSHIP_WNESTING') ###\n" if !$header;
		$header++;
		print STDOUT "$id\n";
	}
}

exit;
