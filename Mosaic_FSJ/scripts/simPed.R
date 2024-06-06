#!/usr/bin/env Rscript

# simPed.R

## parameters
seed <- NA
rp <- 0 # re-paring rate, the rate at members of breeding pairs change not due to death conditional on available breeders (this does not affect the number of nests).
ancfile <- NA # file of input ancestor color IDs
maxoffspring <- 5 # maximum number of fledged offspring observed in Mosaic monitoring data
mature <- 1 # minimum age in years at which indivdiuals can produce offspring
p_male <- 0.5 # probility that an offspring is male for assigning sex to individuals with unknown sex
keep_unbanded = 0 # whether to keep unbanded individuals in analyses
rmatrix = 0 # type of relatedness matrix to calculate (0 = none, 1 = additive)
outprefix = NULL

## parse inputs
parseArgs <- function(argv = NULL) {
	arg <- vector(mode = "list", length=12) # [ped TSV, individual metadata TSV, nest metadata TSV, re-pair rate, ancestral ID file, max number offspring, sexual maturity, keep unbanded, seed]
	arg[[1]] <- argv[length(argv)-2]
	arg[[2]] <- argv[length(argv)-1]
	arg[[3]] <- argv[length(argv)]
	argv  <- argv[-c(length(argv), length(argv)-1, length(argv)-2)]
	i = 1
	while (i < length(argv)) {
		if (argv[i] == "--rp") {
			arg[[4]] <- as.numeric(argv[i+1])
		} else if (argv[i] == "--anc") {
			arg[[5]] <- argv[i+1]
		} else if (argv[i] == "--maxoffspring") {
			arg[[6]] <- as.integer(argv[i+1])
		} else if (argv[i] == "--mature") {
			arg[[7]] <- as.integer(argv[i+1])
		} else if (argv[i] == "--keep_unbanded") {
			arg[[8]] <- 1
			i = i+1
			next
		} else if (argv[i] == "--p_male") {
			arg[[9]] <- as.numeric(argv[i+1])
		} else if (argv[i] == "--rmatrix") {
			arg[[10]] <- as.integer(argv[i+1])
		} else if (argv[i] == "--out") {
			arg[[11]] <- argv[i+1]
		} else if (argv[i] == "--seed") {
			arg[[12]] <- as.integer(argv[i+1])
		} else {
			arg <- paste0("Unknown argument ", argv[i], "\n")
			break
		}
		i = i+2
	}
	return(arg)
}

args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 3) {
	cat("\nsimPed.R [options] <ped tsv file> <individual metadata tsv> <nest metadata>\n")
	cat("\nped TSV file must have a header with the fields:\nID, SIRE_ID, DAM_ID, SEX, COHORT, COHORT_LAST, POPULATION\n")
	cat("\nindividual metadata TSV file must have a header with fields:\nCOLOR_ID, YEAR_CLASS, BIOLOGICAL_ORIGIN_STATUS\n")
	cat("\nnest metadata TSV file must have a header with fields:\nBREEDING_MALE, BREEDING_FEMALE, PROPERTY, YEAR\n")
	cat("\nOptional arguments:\n")
	cat(paste0("--rp <STR|FLOAT>   Rate per generation that breeding pairs are replaced by new pairs, \"e\" empirically reforms pairs based on nest data. [", rp, "]\n"))
	cat("--anc <STRING>   File of ancestor color band IDs (one per row).\n")
	cat(paste0("--maxoffspring <INT>   Max number offspring for a breeding pair per generation. [", maxoffspring, "]\n"))
	cat(paste0("--mature <INT>   Minimum age in years at which an individual can reproduce. [", mature, "]\n"))
	cat("--keep_unbanded   Keeps unbanded individuals in analysis, which are discarded by default.\n")
	cat(paste0("--p_male <FLOAT>   Proportion of individuals of unknown sex that are male. [", p_male, "]\n"))
	cat(paste0("--rmatrix <INT>   Calculate relatedness matrix, 0 = none, 1 = additive. [", rmatrix, "]\n"))
	cat("--out <STRING>   Output file name (prefix). If not specified, output files have '_seed.extension' appended to ped input file name.\n")
	cat("--seed <INT>   Sets the random seed. [random]\n")
	cat("\n")
	options(show.error.messages=FALSE)
	stop()
}

arg <- parseArgs(args)

for (i in 1:3) {
	if (file.exists(arg[[i]][1]) == FALSE) stop(paste0("Unable to find file: ", arg[[i]][1]))
}

ped <- read.table(arg[[1]][1], head=TRUE)
colnames(ped) <- toupper(colnames(ped))
metadat <- read.csv(arg[[2]][1], head=TRUE, sep="\t")
colnames(metadat) <- toupper(colnames(metadat))
nestdat <- read.csv(arg[[3]][1], head=TRUE, sep="\t")
colnames(nestdat) <- toupper(colnames(nestdat))
if (!is.null(arg[4])) rp = arg[[4]][1]
if (is.numeric(rp)) {
	if (rp < 0 || rp > 1) stop("--rp float must be in [0,1]")
} else {
	rp = tolower(rp)
	if (rp != "e") stop(paste0("Invalid --rp ", rp, ", argument can either by 'e' or FLOAT.\n"))
	stop("Empiric re-pairing rates not implemented yet\n") # TL TODO
}
anc = NULL
if (!is.null(arg[[5]])) {
	ancfile = arg[[5]][1]
	anc = read.table(ancfile, head=FALSE)$V1
}
if (!is.null(arg[[6]])) maxoffspring = arg[[6]][1]
if (maxoffspring < 1) stop(paste0("Invalid --maxoffspring ", maxoffspring, ". Max offspring must be a positive integer.\n"))
if (!is.null(arg[[7]])) mature = arg[[7]][1]
if (mature < 0) stop(paste0("Invalid --mature ", mature, ". Age of sexual maturity must be >= 0.\n"))
if (!is.null(arg[[8]])) keep_unbanded = arg[[8]][1]
if (!is.null(arg[[9]])) p_male = arg[[9]][1]
if (p_male > 1 || p_male < 0) stop(paste0("Invalid --p_male ", p_male, ". Value must be in range [0,1].\n"))
if (!is.null(arg[[10]])) rmatrix = arg[[10]][1]

if (rmatrix != 0 & rmatrix != 1) stop(paste0("Invalid --rmatrix ", rmatrix, ". Value can be either 0 or 1\n"))
if (rmatrix == 1) {
	if (!require("nadiv")) stop(paste0("ERROR. --rmatrix ", rmatrix, " requires 'nadiv' package\n"))
}

if (!is.null(arg[length(arg)])) seed = arg[[length(arg)]][1] else seed = round(runif(1,min=1,max=10^9))
set.seed(seed)

if (!is.null(arg[[11]])) outprefix = arg[[11]][1] else outprefix = paste0(arg[[1]][1],"_",seed)


write(paste0("\nSIMPED INPUTS\n--------------\nped file: ", arg[[1]][1], "\nindividual metadata file: ", arg[[2]][1], "\nnest metadata file: ", arg[[3]][1], 
"\nre-pair rate: ", rp, "\nancestral IDs file: ", ancfile, "\nmaximum offspring per pair each generation: ", maxoffspring, 
"\nAge of sexual maturity: ", mature, "\nProportion males: ", p_male, "\nRelatedness matrix: ", rmatrix, "\nseed: ", seed), stderr())
if (keep_unbanded == 1) write("Discarding unbanded individuals\n", stderr()) else write("Keeping unbanded individuals\n", stderr())

## set metadata info

ped$ORIGIN = sapply(ped$ID, function(x,df){df$BIOLOGICAL_ORIGIN_STATUS[which(df$COLOR_ID == x)]}, USE.NAMES=FALSE, df=metadat)

ped$YEAR_CLASS = sapply(ped$ID, function(x,df){df$YEAR_CLASS[which(df$COLOR_ID == x)]}, USE.NAMES=FALSE, df=metadat)
# set year class of "Unknown_Adult" to (COHORT - mature value) to ensure that these individuals can immediately be chosen as breeders
idx = which(ped$YEAR_CLASS == "Unknown_Adult")
ped$YEAR_CLASS[idx] = ped$COHORT[idx] - mature
# set year class of "unknown" to (COHORT - mature value) to ensure that these individuals can immediately be chosen as breeders
# Note that the only such individual in the Mosaic pedigree, "FP-AS", entered the Core Region in 2021 and was banded as AHY in 12/29/2020 
# so treating "unknown" as a sexually mature seems like an okay assumption.
idx = which(ped$YEAR_CLASS == "unknown")
ped$YEAR_CLASS[idx] = ped$COHORT[idx] - mature
ped$YEAR_CLASS <- as.numeric(ped$YEAR_CLASS)
# check that the year class of a local recruit matches its cohort value - if not, force it to be and issue a warning
mismatch_idx = which(ped$ORIGIN == "LR" & ped$COHORT != ped$YEAR_CLASS)
if (length(mismatch_idx) > 0) {
	write(paste0("WARNING. Forcing input YEAR_CLASS to match COHORT for ", paste(ped$ID[mismatch_idx], collapse=", ")), stderr())
	ped$YEAR_CLASS[mismatch_idx] = ped$COHORT[mismatch_idx]
}

ped$FOCAL_ANC <- 0
if (!is.null(anc)) { 
	ped$FOCAL_ANC = sapply(ped$ID, function(x,v){ifelse(x %in% v == TRUE, 1, 0)}, USE.NAMES=FALSE, v=anc)
	# Local recruits denoted as ancestors in anc file are treated as residents, i.e. they will not have parental information 
	# but will contribute to the population size and serve as potential breeders.
	ped$ORIGIN[which(ped$ORIGIN == "LR" & ped$FOCAL_ANC == 1)] = "RESIDENT_LR"
}

# Treat a CR inhabitant born outside of the CR to a CR parent that emigrated out of the CR as "immigrant local recruits".
# Example: HS-BF's offspring. HS-BF was born in the CR, emigrated to Moody, and their offspring F-HSB and F-HRZ immigrated into the CR.
ilr.df = NULL
for (i in 1:nrow(ped)) {
	sire = '*'
	dam = '*'
	if (ped$ORIGIN[i] == "LR" || ped$ORIGIN[i] == "RESIDENT_LR") next
	if (ped$SIRE_ID[i] != '*') {
		sire_idx = which(ped$ID == ped$SIRE_ID[i])
		if (length(sire_idx) > 0) {
			if (ped$ORIGIN[sire_idx] == "LR" || ped$ORIGIN[sire_idx] == "RESIDENT_LR") sire = 1 else sire = 0	
		}
	}
	if (ped$DAM_ID[i] != '*') {
		dam_idx = which(ped$ID == ped$DAM_ID[i])
		if (length(dam_idx) > 0) {
			if (ped$ORIGIN[dam_idx] == "LR" || ped$ORIGIN[dam_idx] == "RESIDENT_LR") dam = 1 else dam = 0		
		}
	}
	if (sire == 1 || dam == 1) {
		ilr.df = rbind(ilr.df, data.frame(ID = ped$ID[i], COHORT = ped$COHORT[i], SIRE_CR = sire, DAM_CR = dam)) # [ID, COHORT, 1=SIRE IS LOCAL RESIDENT, 1=DAM IS LOCAL RESIDENT]
	}
}
ped$ORIGIN[which(ped$ID %in% ilr.df$ID == TRUE)] = "I_LR"

# count types of pairs involving emigrants produced offspring (one or both parents emigrants)
mpair_type_counts <- c(length(which(ilr.df$SIRE_CR == 1 & ilr.df$DAM_CR == 1)), length(which(ilr.df$SIRE_CR == 1 & (ilr.df$DAM_CR == 0 | ilr.df$DAM_CR == '*'))), length(which((ilr.df$SIRE_CR == 0 | ilr.df$SIRE_ID == '*') & ilr.df$DAM_CR == 1)))
migrant_pair_wts = c(mpair_type_counts[1], sum(mpair_type_counts[2:3]))/nrow(ilr.df) # [P(both parents CR emigrants), P(one parent CR emigrant)]

# Check sexes based on whether an individual was a known sire or dam and fix them if necessary
sire_idx <- which(ped$ID %in% ped$SIRE_ID)
badmale <- ped$ID[which(ped$ID %in% ped$SIRE_ID & ped$SEX == "FEMALE")]
if (length(badmale) > 0) write(paste0("WARNING. Found female sires in ped input (forcing them to be male): ", paste(badmale, collapse=", ")), stderr())
ped$SEX[sire_idx] = "MALE"
dam_idx <- which(ped$ID %in% ped$DAM_ID)
badfemale <- ped$ID[which(ped$ID %in% ped$DAM_ID & ped$SEX == "MALE")]
if (length(badfemale) > 0) write(paste0("WARNING. Found male dams in ped input (forcing them to be female): ", paste(badfemale, collapse=", ")), stderr())
ped$SEX[dam_idx] = "FEMALE"

# Check for individuals that had LR offspring after they were dead (after pedigree COHORT_LAST which was the last time 
# they were observed in the CR). These individuals likely were alive in the CR but missed in the census so set their 
# COHORT_LAST to their latest LR COHORT. 
# Also log information for emigration events by CR inhabitants that contributed offspring to the CR from outside of it (offspring 
# need to be immigrants into the CR).

ped$EMIGRATE = 0

adjSurvival <- function(ped = NULL, pcol = NA) {
	# ped: pedigree dataframe
	# pcol: column of parent in ped dataframe
	for (i in 1:nrow(ped)) {
		maxyc_i = 0
		maxyc_lr = 0
		p_idx = NULL
		if (ped[i,pcol] != "*") p_idx = which(ped$ID == ped[i,pcol])
		if (length(p_idx) > 0) {
			i_offspring_idx = which(ped[,pcol] == ped[i,pcol] & ped$ORIGIN == "I_LR")
			if (length(i_offspring_idx) > 0) maxyc_i = max(ped$YEAR_CLASS[i_offspring_idx])
			lr_offspring_idx = which(ped[,pcol] == ped[i,pcol] & (ped$ORIGIN %in% c("LR", "RESIDENT_LR") == TRUE ))
			if (length(lr_offspring_idx) > 0) maxyc_lr = max(ped$YEAR_CLASS[lr_offspring_idx])
			if (maxyc_lr > ped$COHORT_LAST[p_idx]) {
				write(paste0("WARNING. ",ped$ID[p_idx], " produced LR offspring after last seen in census. Adjusting COHORT_LAST."), stderr())
				ped$COHORT_LAST[p_idx] = maxyc_lr
			}
			if (maxyc_i > ped$COHORT_LAST[p_idx]) {
				ped$EMIGRATE[p_idx] = ped$COHORT_LAST[p_idx] + 1 # This is the year an individual migrates out of the CR
				ped$COHORT_LAST[p_idx] = maxyc_i # COHORT_LAST is now the last known year of life for individual outside of the CR
			}
		}
	}
	return(ped)
}

ped <- adjSurvival(ped=ped, pcol=which(colnames(ped) == "SIRE_ID"))
ped <- adjSurvival(ped=ped, pcol=which(colnames(ped) == "DAM_ID"))

# set probability that an emigrant will be male or female
total_emigrants = length(which(ped$EMIGRATE != 0))
p_male_migrant = ifelse(total_emigrants > 0,length(which(ped$SEX == "MALE" & ped$EMIGRATE != 0))/total_emigrants, 0.5) # probability that CR emigrant is male

# Assign sex to individuals with unknown sex
# randomly assign sex according to input probability of being male
mis_sex = which(ped$SEX != "MALE" & ped$SEX != "FEMALE")
ped$SEX[as.numeric(sample(as.character(mis_sex), size=as.integer(length(mis_sex)*p_male), replace=FALSE))] = "MALE"
ped$SEX[which(ped$SEX != "MALE" & ped$SEX != "FEMALE")] = "FEMALE"

# limit breeding pairs to those in the Core Region
nestdat <- nestdat[which(nestdat$PROPERTY == "MW" | nestdat$PROPERTY == "COKER" | nestdat$PROPERTY == "DUETTE" | nestdat$PROPERTY == "E COKER A"),]
nestdat$PAIR <- paste(nestdat$BREEDING_MALE, nestdat$BREEDING_FEMALE, sep="_")

if (keep_unbanded == 0) {
	# remove unbanded info from pedigree     
	idx = grep("unband|unb|none|noband|unknown", ped$ID, ignore.case=TRUE)
	if (length(idx) > 0) ped <- ped[-idx,]
	# set ubanded ancestors to missing
	idx = grep("unband|unb|none|noband|unknown", ped$SIRE_ID, ignore.case=TRUE)
	if (length(idx) > 0) ped$SIRE_ID[idx] <- '*'                                            
	idx = grep("unband|unb|none|noband|unknown", ped$DAM_ID, ignore.case=TRUE)
	if (length(idx) > 0) ped$DAM_ID[idx] <- '*'

	# remove unbanded individuals from indivdiual metadata
	idx = grep("unband|unb|none|noband|unknown", metadat$COLOR_ID, ignore.case=TRUE)
	if (length(idx) > 0) metadat <- metadat[-idx,]

	# remove nests involving unbanded individuals
	idx = sort(unique(c(grep("unband|unb|none|noband|unknown", nestdat$BREEDING_MALE, ignore.case=TRUE), grep("unband|unb|none|noband|unknown", nestdat$BREEDING_FEMALE, ignore.case=TRUE))))
	if (length(idx) > 0) nestdat <- nestdat[-idx,]
}

## calculate empirical parameters
# estimate empirical mate switching rates (not due to mate death) here. TL TODO

## simulate pedigree

ped.sim = NULL
tmin = min(ped$COHORT)
tmax = max(ped$COHORT)
breeders = NULL # all potential breeders in the focal population
bp = data.frame(matrix(NA,nrow=length(unique(nestdat$PAIR)), ncol=3)) # all breeding pairs [male, female, number offspring in generation]
colnames(bp) = c("male", "female", "noffspring")
bp[,1] = as.character(bp[,1])
bp[,2] = as.character(bp[,2])
bp[,3] = as.numeric(bp[,3])
bp_migrant = NULL
anc_id = NULL # IDs of ancestral individuals in simulated pedigree
id_map = matrix(NA, nrow=length(unique(c(ped$ID, ped$SIRE_ID, ped$DAM_ID))) - length(which(ped$ORIGIN == "LR")), ncol=2) # maps original color band IDs to ID in simulation for non-local recruits (i.e. residents, translocated individuals, immigrants)
alive_migrant = NULL
map_c = 1
id = 1

# Assign sim IDs to all parents not part of the core population - could use this to assign family structure to ancestors with unknown relatedness.
# Will keep unknown parents of founders as '*', i.e. unknown.
nc_parents = unique(c(ped$SIRE_ID[which(ped$SIRE_ID %in% ped$ID == FALSE)], ped$DAM_ID[which(ped$DAM_ID %in% ped$ID == FALSE)]))
nc_parents = nc_parents[-which(nc_parents == '*')]
for (pid in nc_parents) {
	id_map[map_c,] = c(pid, id)
	map_c = map_c+1
	id = id+1
}

for (t in tmin:tmax) {
	# Insert non-local recruit individuals into the population.
	# RESIDENT_LR are treated as immigrants that spontaneously appear in the population while maintaining parent-offspring relationships to other founders
	# note: pedidx defined suched that RESIDENT_LR individuals are added last, which ensures that their parents will be in the id_map
	pedidx = c(which(ped$COHORT == t & ped$ORIGIN != "LR" & ped$ORIGIN != "I_LR" & ped$ORIGIN != "RESIDENT_LR"), which(ped$COHORT == t & ped$ORIGIN == "RESIDENT_LR"))
	addstart = ifelse(is.null(ped.sim), 1, nrow(ped.sim)+1)
	ped.sim = rbind(ped.sim, ped[pedidx,])
	for (i in addstart:nrow(ped.sim)) {
		# assign sim ID to focal individual
		map_idx = which(id_map[,1] == ped.sim$ID[i])
		if (length(map_idx) > 0) {
			ped.sim$ID[i] = id_map[map_idx,2]
		} else {
			id_map[map_c,] = c(ped.sim$ID[i], id)
			map_c = map_c+1
			ped.sim$ID[i] = id
			id = id+1
		}

		# Set parent info of new individual. Keeping missing parent info as '*'.
		if (ped.sim$SIRE_ID[i] != '*') {
			map_idx = which(id_map[,1] == ped.sim$SIRE_ID[i])
			if (length(map_idx) > 0) {
				ped.sim$SIRE_ID[i] = id_map[map_idx,2]
			} else {
				id_map[map_c,] = c(ped.sim$SIRE_ID[i], id)
				map_c = map_c+1
				ped.sim$SIRE_ID[i] = id
				id = id+1
			}
		}
		if (ped.sim$DAM_ID[i] != '*') {
			map_idx = which(id_map[,1] == ped.sim$DAM_ID[i])
			if (length(map_idx) > 0) {
				ped.sim$DAM_ID[i] = id_map[map_idx,2]
			} else {
				id_map[map_c,] = c(ped.sim$DAM_ID[i], id)
				map_c = map_c+1
				ped.sim$DAM_ID[i] = id
				id = id+1
			}
		}
	}

	# Manage survival
	# update individuals alive in the population
	alive <- ped.sim[which(ped.sim$COHORT_LAST >= t),] # all individuals alive in the CR
	if (!is.null(alive_migrant) && nrow(alive_migrant) > 0) {
		midx = which(alive$ID %in% alive_migrant$ID)
		if (length(midx) > 0) alive = alive[-midx,]	
	}
	if (nrow(alive) < 1) {
		write(paste0("WARNING. Population is extinct at time", t), stderr())
		break
	}

	# Manage migration
	if (!is.null(alive_migrant)) {
		rmv_midx = which(alive_migrant$COHORT_LAST < t)
		if (length(rmv_midx) > 0) alive_migrant = alive_migrant[-rmv_midx,]
	}
	mig_n = length(which(alive$EMIGRATE == t))
	if (mig_n > 0) {
		mig_avail = which((t - alive$YEAR_CLASS) >= 1 & alive$ID %in% c(bp[,1],bp[,2]) == FALSE ) # individuals need to be at least 1 year old and not in a breeding pair to migrate
		nemigrate = 0
		if (length(mig_avail) > 0) {
			wt = alive$SEX[mig_avail]
			wt = replace(wt,which(wt == "MALE"),p_male_migrant)
			wt = replace(wt,which(wt == "FEMALE"), 1-p_male_migrant)
			wt = as.numeric(wt)
			nemigrate = min(mig_n, length(which(wt > 0)))
			if (nemigrate > 0) {
				migrants = as.numeric(unique(sample(as.character(mig_avail), size=nemigrate, replace=FALSE, prob=wt)))
				alive_migrant <- rbind(alive_migrant, alive[migrants,])
				ped.sim[which(ped$ID == alive$ID[migrants]), which(colnames(ped) %in% c("COHORT_LAST", "EMIGRATE") == TRUE)] = t 
				alive <- alive[-migrants,]
			}
		}
		# if there are not enough qualifying migrants, randomly kill off individuals from the population to match empirical population size
		nremove = n_mig - nemigrate
		if (nremove > 0) {
			nremove = min(nremove, nrow(alive))
			rmv_idx = as.numeric(unique(sample(as.character(1:nrow(alive)), size=nremove, replace=FALSE)))
			ped.sim$COHORT_LAST[which(ped$ID == alive$ID[rmv_idx])] = t
			alive = alive[-rmv_idx,]
		}
	}

	# update potential breeders
	breeders.m <- alive$ID[which((t - alive$YEAR_CLASS) >= mature & alive$SEX == "MALE")] # males capable of reproducing
	breeders.f <- alive$ID[which((t - alive$YEAR_CLASS) >= mature & alive$SEX == "FEMALE")] # females capabale of reproducing

	# Manage breeders
	# first, dissipate breeding pairs containing a dead individual.
	rmvidx <- which((bp[,1] %in% alive$ID == FALSE | bp[,2] %in% alive$ID == FALSE) & !is.na(bp[,1]))
	if (length(rmvidx) > 0) bp[rmvidx,] <- NA

	# mask individuals available for breeding if they are already in a pair
	breeders.m.avail = breeders.m
	if (length(breeders.m.avail) > 0) {
		midx = which(breeders.m.avail %in% bp[,1])
		if (length(midx) > 0) breeders.m.avail <- breeders.m.avail[-midx]
	}
	breeders.f.avail = breeders.f
	if (length(breeders.f.avail) > 0) {
		fidx = which(breeders.f.avail %in% bp[,2])
		if (length(fidx) > 0) breeders.f.avail <- breeders.f.avail[-fidx]
	}
	
	# Adjust number of breeding pairs to match current generation empirical number
	nests.n = length(unique(nestdat$PAIR[which(nestdat$YEAR == t)])) # number empirical pairs in existence
	
	#number of empirical pairs that only allows an individual to appear in one pair per season
	#nestsub = nestdat[which(nestdat$YEAR == t),]
	#subidx = sapply(unique(nestsub$BREEDING_MALE), function(x,df){which(df$BREEDING_MALE == x)[1]}, df=nestsub, USE.NAMES=FALSE)
	#nests.n = length(unique(nestsub$BREEDING_FEMALE[subidx]))	

	simpair.n = length(which(!is.na(bp[,1]))) # number of simulated pairs in existence

	# add breeding pairs as necessary as long as there are available breeders
	while (simpair.n < nests.n && length(breeders.m.avail) > 0 && length(breeders.f.avail) > 0) {
		midx = sample(1:length(breeders.m.avail), 1, replace=FALSE)
		fidx = sample(1:length(breeders.f.avail), 1, replace=FALSE)
		emptyidx = which(is.na(bp[,1]))
		if (length(emptyidx) > 0) {
			bp[emptyidx[1],] = list(breeders.m.avail[midx], breeders.f.avail[fidx], NA)
		} else {
			# expand breeding pair matrix if necessary
			bp[nrow(bp)+1,] = NA
			bp[nrow(bp),] = list(breeders.m.avail[midx], breeders.f.avail[fidx], NA)
		}
		breeders.m.avail <- breeders.m.avail[-midx]
		breeders.f.avail <- breeders.f.avail[-fidx]
		simpair.n = simpair.n+1
	}

	# remove breeding pairs as necessary
	while (simpair.n > nests.n) {
		delidx = sample(which(!is.na(bp[,1])), 1, replace=FALSE)
		breeders.m.avail = c(breeders.m.avail, bp[delidx,1])
		breeders.f.avail = c(breeders.f.avail, bp[delidx,2])
		bp[delidx,] <- NA
		simpair.n = simpair.n-1
	}

	# replace individuals comprising breeding pairs according to repairing rate
	pairidx = which(!is.na(bp[,1]))
	simpair.n = length(pairidx)
	if (simpair.n > 0 && rp > 0) {
		i = 1
		rv = runif(2*simpair.n)
		for (idx in pairidx) {
			if (rv[i] <= rp && length(breeders.m.avail) > 0) {
				midx = sample(1:length(breeders.m.avail), 1)
				m.id = bp[idx,1]
				bp[idx,1] = breeders.m.avail[midx]
				breeders.m.avail[midx] = m.id
			}
			i = i+1
			if (rv[i] <= rp && length(breeders.m.avil) > 0) {
				fidx = sample(1:length(breeders.f.avail), 1)
				f.id = bp[idx,2]
				bp[idx,2] = breeders.f.avail[fidx]
				breeders.f.avail[fidx] = f.id
			}
		}
	}

	# Update pool of emigrant breeders. This consists of pairs where at least 1 individual emigrated from the CR.
	# Add individuals back to pool of available migrant breeders if their mate is dead.
	if (!is.null(bp_migrant) && nrow(bp_migrant) > 0) {
		badpair_idx = which((bp_migrant[,1] != '*' & bp_migrant[,1] %in% alive_migrant == FALSE) | (bp_migrant[,2] != '*' & bp_migrant[,2] %in% alive_migrant == FALSE))
		if (length(badpair_idx) > 0) bp_migrant = bp_migrant[-badpair_idx,]
	}
	# form new migrant breeding pairs
	if (!is.null(alive_migrant) && nrow(alive_migrant) > 0) {
		avail_idx = which((t - alive_migrant$YEAR_CLASS) >= mature & alive_migrant$ID %in% c(bp_migrant[,1], bp_migrant[,2]) == FALSE)
		m_avail = alive_migrant[avail_idx,]
		if (length(avail_idx) > 0) {
			# sample having 2 or 1 CR individuals in a pair. Sample number is max number of pairs if each pair has only 1 CR individual.
			ptypes = sample(c(2,1), size=nrow(m_avail), replace=TRUE, prob=migrant_pair_wts)
			for (k in ptypes) {
				if (k == 2) {
					midx = which(m_avail$SEX == "MALE")
					fidx = which(m_avail$SEX == "FEMALE")
					if (length(midx) > 0 && length(fidx) > 0) {
						msamp = as.numeric(sample(as.character(midx), size=1, replace=FALSE))
						fsamp = as.numeric(sample(as.character(fidx), size=1, replace=FALSE))
						bp_migrant = rbind(bp_migrant, data.frame(male=m_avail$ID[msamp], female=m_avail$ID[fsamp]))
						m_avail = m_avail[-c(msamp, fsamp)]
					}
				} else {
					isamp = as.numeric(sample(as.character(1:nrow(m_avail)), size=1, replace=FALSE))
					if (m_avail$SEX[isamp] == "MALE") {
						bp_migrant = rbind(bp_migrant, data.frame(male=m_avail$ID[isamp], female='*', noffspring = as.numeric(NA)))
					} else bp_migrant = rbind(bp_migrant, data.frame(male='*', female=m_avail$ID[isamp], noffspring = as.numeric(NA)))
					m_avail = m_avail[-isamp,]
				}
				if (nrow(m_avail) < 1) break
			}
		}
	}

	# Manage reproduction
	# Number of individuals born into the population. Individuals treated as founding ancestors in the anc file are excluded from this count 
	# including local recruits denoted as founding ancestors.
	recruit_idx = NULL
	lr_idx = which(ped$COHORT == t & ped$ORIGIN == "LR")
	ilr_idx = which(ped$COHORT == t & ped$ORIGIN == "I_LR")
	all_lr_idx = c(lr_idx, ilr_idx)
	offspring_total = length(lr_idx) # number of LR offspring to assign parents to in current generation
	bpidx = which(!is.na(bp[,1]))
	if (length(bpidx) > 0 && offspring_total > 0) {
		recruit_idx = as.numeric(sample(as.character(all_lr_idx), size=offspring_total, replace=FALSE)) # ensures random assignment of LR and I_LR offsprings to parents
		bp[bpidx,3] = 0
		#bpvec = sample(x = as.character(bpidx), size = offspring_total, replace = TRUE) # as.character to prevent 1:x, also consider adding weights according to parent age		
		nassign = offspring_total - sum(bp[,3], na.rm=TRUE) # number of offspring needing to be assigned a breeding pair

		while (nassign > 0 && length(bpidx) > 0) {
			nassign = offspring_total - sum(bp[,3], na.rm=TRUE)
			bpvec = sample(x = as.character(bpidx), size = nassign, replace = TRUE) # as.character to prevent 1:x, also consider adding weights according to parent age
			offspring_dist = table(bpvec)
			bp[as.numeric(names(offspring_dist)),3] = bp[as.numeric(names(offspring_dist)),3] + unname(offspring_dist)
			bp[which(bp[,3] > maxoffspring),3] = maxoffspring # limit maximum offspring
			bpidx = which(!is.na(bp[,1]) & bp[,3] < maxoffspring) # remaining parents that can have more offspring
			nassign = offspring_total - sum(bp[,3], na.rm=TRUE)
		}
	}

	# Update pedigree with LR individuals
	bp_success_idx = which(bp[,3] > 0)
	n = 0
	if (length(bp_success_idx) > 0) n = sum(bp[bp_success_idx,3])
	if (n > 0) {
		addstart = ifelse(nrow(ped.sim) > 0, nrow(ped.sim)+1, 1)
		ped.sim = rbind(ped.sim, ped[recruit_idx,])
		ped.sim$ID[addstart:nrow(ped.sim)] = id:(id+n-1) # assign offspring IDs
		id = id + n
		ped.sim$SIRE_ID[addstart:nrow(ped.sim)] = rep(bp[bp_success_idx,1],bp[bp_success_idx,3]) # assign father ID
		ped.sim$DAM_ID[addstart:nrow(ped.sim)] = rep(bp[bp_success_idx,2],bp[bp_success_idx,3]) # assign mother ID
		ped.sim$ORIGIN[addstart:nrow(ped.sim)] = "LR"
	}

	# insert I_LR individuals into pedigree. These are individuals born to CR emigrants that imigrate into the CR.
	recruit_idx = all_lr_idx[which(all_lr_idx %in% recruit_idx == FALSE)]
	n_mig = length(recruit_idx)
	if (!is.null(bp_migrant) && nrow(bp_migrant) > 0) bp_migrant[,3] = 0
	for (i in recruit_idx) {
	# adding migrants like this is perhaps a bit inefficient, but this is okay because don't expect many I_LR individuals per generation
		sire_id = '*'
		dam_id = '*'
		if (!is.null(bp_migrant) && nrow(bp_migrant) > 0) {
			bpm_idx = which(bp_migrant[,3] < maxoffspring)
			if (length(bpm_idx) > 0) {
				mpair_idx = as.numeric(sample(as.character(bpm_idx), size=1, replace=FALSE))
				sire_id = bp_migrant[mpair_idx,1]
				dam_id = bp_migrant[mpair_idx,2]
				bp_migrant[mpair_idx,3] = bp_migrant[mpair_idx,3]+1
			}
		}
		ped.sim = rbind(ped.sim, ped[i,])
		ped.sim$ID[nrow(ped.sim)] = id
		id = id+1
		ped.sim$SIRE_ID[nrow(ped.sim)] = sire_id
		ped.sim$DAM_ID[nrow(ped.sim)] = dam_id
		ped.sim$ORIGIN[nrow(ped.sim)] = "I_LR"
	}

}

# Set COHORT_LAST for input pedigree emigrants to year before emigration since this was the last time they were in the CR
ped_emig_idx = which(ped.sim$EMIGRATE != 0)
ped.sim$COHORT_LAST[ped_emig_idx] = ped.sim$EMIGRATE[ped_emig_idx]-1

# Recode population IDs for clearer downstream tracking
ped.sim$POPULATION[which(ped.sim$ORIGIN == "IE" | ped.sim$ORIGIN == "RESIDENT_LR")] = "CORE_RESIDENT"
ped.sim$POPULATION[which(ped.sim$ORIGIN == "LR")] = "CORE_LR"
ped.sim$POPULATION[which(ped.sim$ORIGIN == "I_LR")] = "UNKNOWN"

#print(ped.sim); stop("Not actually an error, debug stop point\n") # debug

# write pedigree file
outped = paste0(outprefix,".ped")
write.table(ped.sim[,c("ID","SIRE_ID","DAM_ID","SEX","COHORT","COHORT_LAST","POPULATION")], file=outped, col.names=TRUE, row.names=FALSE, sep="\t", quote=FALSE)
write(paste0("Output simulated pedigree: ",outped), stderr())

# make dataframe of focal individual ID map
if (!is.null(anc)) {
	id_map_focal = data.frame(SIM_ID = ped.sim$ID[which(ped.sim$FOCAL_ANC == 1)])
	id_map_focal$PED_ID = sapply(id_map_focal$SIM_ID, function(x,df){df[which(df[,2] == x),1]}, df=id_map, USE.NAMES=FALSE)
	out_idmap = paste0(outprefix, ".id")
	write.table(id_map_focal, file=out_idmap, col.names=TRUE, row.names=FALSE, sep="\t", quote=FALSE)
	write(paste0("Output ancestral individual ID map: ", out_idmap), stderr())
}

# calculate additive relatedness matrix
if (rmatrix == 1) {
	ped2 = prepPed(ped.sim[,1:4], gender='SEX', check=TRUE)
	rmat = as.matrix(makeA(ped2[,1:3])) # as.matrix to convert sparse dsCMatrix into base matrix class
	outmat = paste0(outprefix,".mat")
	write.table(rmat, file=outmat, col.names=TRUE, row.names=FASE, quote=FALSE, sep=" ")
	write(paste0("Output additive relatedness matrix: ",rmat), stderr())
}
