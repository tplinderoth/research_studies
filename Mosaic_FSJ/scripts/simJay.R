#!/usr/bin/env Rscript

# simJay.R

## functions

parseArgs <- function(argv = NULL, params) {
	args = params
	err = 0
	i = 1
	while (i <= length(argv)) {
		if (err > 0) break
		if (argv[i] == "--ind_file") {
			if (length(argv) >= i+1) args$ind_file = argv[i+1] else err = 1
			if (!file.exists(args$ind_file)) {
				write(paste0("Error: ", args$ind_file, " does not exist"), stderr())
				err = 1
			}
		} else if (argv[i] == "--pair_param") {
			if (length(argv) >= i+2) args$pair_param = as.numeric(argv[c(i+1,i+2)]) else err = 1
			i = i+3
			next
		} else if (argv[i] == "--offspring_file") {
			if (length(argv) >= i+1) args$offspring_file = argv[i+1] else err = 1
			if (!file.exists(args$offspring_file)) {
				write(paste0("Error: ", args$offspring_file, " does not exist"), stderr())
				err = 1
			}
		} else if (argv[i] == "--survive_file") {
			if (length(argv) >= i+1) args$survive_file = argv[i+1] else err = 1
			if (!file.exists(args$survive_file)) {
				write(paste0("Error: ", args$survive_file, " does not exist"),stderr())
				err = 1
			}
		} else if (argv[i] == "--male_p") {
			if (length(argv) >= i+1) args$male_p = as.numeric(argv[i+1]) else err = 1
			if (args$male_p < 0 | args$male_p > 1) {
				write(paste0("Error: male_p ", args$male_p, " is out of range [0,1]"),stderr())
				err = 1
			}
		} else if (argv[i] == "--mature") {
			if (length(argv) >= i+1) args$mature = as.integer(argv[i+1]) else err = 1
		} else if (argv[i] == "--max_offspring") {
			if (length(argv) >= i+1) args$max_offspring = as.integer(argv[i+1]) else err = 1
			if (args$max_offspring < 1) {
				write(paste0("Error: --max_offspring", args$max_offspring, " is less than 1"),stderr())
				err = 1
			}
		} else if (argv[i] == "--seed") {
			if (length(argv) >= i+1) args$seed = as.integer(argv[i+1]) else err = 1
		} else if (argv[i] == "--unif_reproduction") {
			args$unif_reproduction = TRUE
			i = i+1
			next
		} else if (argv[i] == "--n_timesteps") {
			if (length(argv) >= i+1) args$n_timesteps = as.integer(argv[i+1]) else err = 1
			if (args$n_timesteps < 1) {
				write(paste0("Error: Invalid number of time steps, ", args$n_timesteps),stderr())
				err = 1
			}
		} 
		else if (argv[i] == "--empiric_survive") {
			args$empiric_survive = TRUE
			i = i+1
			next
		} else if (argv[i] == "--min_survive_obs") {
			if (length(argv) >= i+1) args$min_survive_obs = as.integer(argv[i+1]) else err = 1
			if (args$min_survive_obs < 1){
				write(paste0("Error: --min_survive_obs must be an integer greater than zero\n"), stderr())
				err = 1
			}
		} else if (argv[i] == "--out") {
			if (length(argv) >= i+1) args$out = argv[i+1] else err = 1
		} else if (argv[i] == "--events_file") {
			if (length(argv) >= i+1) args$events_file = argv[i+1] else err = 1
			if (!file.exists(args$events_file)) {
				write(paste0("Error: ", args$events_file, " does not exist"), stderr())
				err = 1
			}
		} else {
			write(paste0("Unknown argument ", argv[i]),stderr())
			err = 1
		}
		i = i+2
	}
	if (err > 0) args = NULL
	return (args)
}

fmtProbDist <- function(dat) {
	# sets probability bounds for categories for easy event calculations
	# take dataframe with columns (1) event (2) probability
	dat[,2] = dat[,2]/sum(dat[,2]) # ensure probabilities sum to 1
	dat$a = c(0, rep(NA,nrow(dat)-1))
	dat$b = c(dat[1,2], rep(NA,nrow(dat)-1))
	for (i in 2:nrow(dat)) {
		dat$a[i] = dat$b[i-1]              
		dat$b[i] = dat$a[i] + dat[i,2]
	}
	return(dat)
}

## parameters

defaults = list()

defaults$seed <- NA # random number seed
defaults$male_p = 0.5 # probability that individual is male
defaults$offspring_file = NA # TSV file with columns (1) Number fledglings (2) probability
defaults$pair_param = c(0,0) # space-delimited intercept and effect of available breeders for poisson model of new breeding pairs as a function of available breeders (beta0 beta_breeders)
defaults$survive_file = NA # TSV file with columns (1) Number years alive (2) probability
defaults$empiric_survive = FALSE
defaults$min_survive_obs = 1 # Minimum number of survival observations to use empirical survival distribution for a season, otherwise switch to survive_file if it's provided
defaults$ind_file = NA # TSV file with columns ID, SEX, YEAR_CLASS, COHORT, COHORT_LAST, ORIGIN, ANCESTOR, SIRE_ID (optional), DAM_ID (optional)
defaults$mature = 2 # Minimum age at which individual can reproduce
defaults$unif_reproduction = FALSE # Redistribute offspring uniformly randomly across breeding pairs
defaults$n_timesteps = NA # Number of time steps (generations) to simulate
defaults$max_offspring = 5 # maximum number of offspring a breeding pair can have in one season
defaults$out = NA # output file prefix
defaults$events_file = NA # file specifying events

## main

## parse input

args = commandArgs(trailingOnly=TRUE)

if (length(args) < 1) {
	cat ("\nsimJay.R <arguments>\n")
	cat("\n*** Required input ***\n\n")
	cat("--ind_file  <FILE>  TSV file with header and fields: ID, SEX, YEAR_CLASS, COHORT, COHORT_LAST, ORIGIN, ANCESTOR STATUS, SIRE_ID (optional), DAM_ID (optional)\n")
	cat("  ID: individual ID\n  SEX: Male or Female\n  YEAR_CLASS: Year that individual is born\n  COHORT: Year they enter the population\n")
	cat("  COHORT_LAST: Year of death or emmigration from the population\n  ORIGIN: IE, T, LR, I, TD-FR\n  ANCESTOR STATUS: 1 = individual is a population founder, 0 = not a founder\n")
	cat("  SIRE_ID: father identity\n  DAM_ID: mother identity\n")
	cat("  '*' is treated as missing data\n\n")
	cat("--offspring_file  <FILE>  TSV file with fields (1) Number fledglings per pair in a season (2) probability. Lines preceeded with # are ignored\n")
	cat("--pair_param  <FLOAT FLOAT>  Space-delimited parameter values for Poisson model of pair formation: (1) intercept, (2) effect of available breeders. [", defaults$pair_param, "]\n")
	cat("--n_timesteps  <INT>  Number of time steps (generations) to simulate\n")
	cat("--out  <STRING>  Output file prefix\n")
	cat("\n*** Optional input ***\n\n")
	cat("--mature  <INT>  Minimum age at which individuals are eligible to pair and reproduce [", defaults$mature, "]\n")
	cat("--max_offspring  <INT>  Maximum number of offspring that a breeding pair can have in one season [", defaults$max_offspring, "]\n")
	cat("--survive_file  <FILE>  TSV file with fields (1) Number of years alive (2) probability. Lines preceeded with # are ignored\n")
	cat("--empiric_survive  <>  If specified, empirical survival probabilities identified for non-ancestor LRs from ind_file when possible.\n")
	cat("--min_survive_obs  <INT>  Minimum number of survival observations required for empiric distribution, otherwise switch to survive_file distribution if provided [",
	defaults$min_survive_obs, "]\n")
	cat("--events_file  <FILE>  TSV file with head and columns [1] TIME (1-based), [2, ..., N] event parameters \n")
	cat("--male_p  <FLOAT>  Controls the sex ratio and is the probability that an unsexed individual is male [", defaults$male_p, "]\n")
	cat("--unif_reproduction  <>  Evoke to redistribute offspring uniformly at random among breeding pairs during reproduction\n")
	cat("--seed  <INT>  Random number seed [random]\n")
	cat("\n ** Events (specified by header names) **\n")
	cat("TIME:  <INT>  1-based time step at which events happens\n")
	cat("PAIR_LIMIT:  <INT>   No pairing allowed if number pairs exceeds INT ('NA' or '*' values sets no restrictions)\n")
	cat("\n** Input Notes **\n")
	cat("\n* If --empiric_survive is used but the survival distribution cannot be calculated for a time step, survival is assigned based on the survive_file.\n  If survival cannot be empirically calculated and there is no survive file, an error is thrown.\n")
	cat("\n* Without --empiric_survive specified, a --survive_file is required.\n")
	cat("\n* The default --pair_parm, (0,0), specifies no formation of new pairs.\n")
	cat("\n")
	options(show.error.messages=FALSE)
	stop()
}

params = parseArgs(args, defaults)

## check for required inputs

if (is.na(params$ind_file)) stop("Must supply --ind_file\n")
if (is.na(params$offspring_file)) stop("Must supply --offspring_file\n")
if (is.na(params$n_timesteps)) stop("Must specify --n_timesteps\n")
if (params$empiric_survive == FALSE && is.na(params$survive_file)) stop("Must supply survive_file if empiric_survive is false\n");
if (is.na(params$out)) stop("Must supply outfile prefix name with --out\n")

## define outfile names

out_ped = paste0(params$out,".ped") # pedigree output file
out_idmap = paste0(params$out, ".id") # ID map output file
out_stats = paste0(params$out, ".stats") # population stats (e.g. size, number pairs) output file

## set seed

if (is.na(params$seed)) params$seed = round(runif(1,min=1,max=10^9))
set.seed(params$seed)

## adjust parameters if necessary

if (params$empiric_survive == TRUE && is.na(params$survive_file) && params$min_survive_obs > 1) {
	write("--empiric_survive set without --survive_file ==> setting --min_survive_obs to 1", stderr())
	params$min_survive_obs = 1
}

## print parameters for run

write(paste0("\nSIMJAY INPUTS\n--------------\nind_file: ", params$ind_file, "\noffspring probability file: ", params$offspring_file, 
"\npairing parameters: beta0 = ", params$pair_param[1], ", beta1 = ", params$pair_param[2], "\nnumber of timesteps: ", params$n_timesteps, 
"\nsurvival probability file: ", params$survive_file, "\nuse empirical survival probabilities from ind_file: ", params$empiric_survive, 
"\nMinimum number of observations for empirical survival calculation: ", params$min_survive_obs,"\nage of sexual maturity: ", params$mature, 
"\nMax offspring for a breeding pair in a season: ", params$max_offspring, "\nEvents parameter file: ", params$events_file,
"\nproportion males: ", params$male_p, "\ndistribute offspring uniformly at random among breeding pairs: ", params$unif_reproduction, "\nseed: ", params$seed), stderr())

## read input files and set up dataframes

# dataframe individuals to seed population with (and derive empirical survival from)

ind.df <- read.table(params$ind_file, head=TRUE, sep="\t")
colnames(ind.df) = sapply(colnames(ind.df), toupper, USE.NAMES=FALSE)
ind.df = ind.df[order(ind.df$COHORT),]
ind.df$ORIGIN = toupper(ind.df$ORIGIN)
ind.df$YEAR_CLASS = toupper(ind.df$YEAR_CLASS)

# format and set sex
ind.df$SEX = sapply(ind.df$SEX, toupper, USE.NAMES=FALSE)
mis.idx = which(ind.df$SEX == "*")
if (length(mis.idx) > 0) {
	ind.df$SEX[mis.idx] = sample(c('MALE','FEMALE'),size=length(mis.idx), replace=TRUE, prob=c(params$male_p, 1-params$male_p))
}

# format and set age information
ind.df$YEAR_CLASS[which(ind.df$YEAR_CLASS == "*")] = NA
if (length(which(is.na(ind.df$COHORT))) > 0) stop("There is missing COHORT information in ind_file\n")
if (length(which(is.na(ind.df$COHORT_LAST))) > 0) stop("There is missing COHORT_LAST information in ind_file\n")
adult_idx = which(ind.df$YEAR_CLASS == "ADULT")
ind.df$YEAR_CLASS[adult_idx] = ind.df$COHORT[adult_idx] - params$mature # this hack means that all incoming adults of unknown age will be eligible to pair and breed immediately upon entering the population
immi_idx = which((is.na(ind.df$YEAR_CLASS) | ind.df$YEAR_CLASS == "UNKNOWN") & (ind.df$ORIGIN == "I" | ind.df$ORIGIN == "TD-FR"))
ind.df$YEAR_CLASS[immi_idx] = ind.df$COHORT[immi_idx] - params$mature # assume that immigrants with unknown age are mature adults capable of breeding (does not apply to translocated individuals)
ind.df$YEAR_CLASS[which(ind.df$YEAR_CLASS == "UNKNOWN")] = NA
na.idx = which(is.na(ind.df$YEAR_CLASS))
if (length(na.idx) > 0) {
	mis_age = ind.df$ID[which(is.na(ind.df$YEAR_CLASS) & ind.df$ANCESTOR == 1)]
	write(paste0("Warning. Age class information for the following ancestors is missing and so will be set to their input COHORT value: ", mis_age))
	ind.df$YEAR_CLASS[na.idx] = ind.df$COHORT[na.idx] # assume all remaining individuals with unknown age are immature
}
ind.df$YEAR_CLASS = as.numeric(ind.df$YEAR_CLASS)

# format origin info
ind.df$ORIGIN[which(ind.df$ORIGIN == "T")] = "I" # treat translocated individuals as immigrants
ind.df$ORIGIN[which(ind.df$ORIGIN == "TD-FR")] = "I"
ind.df$ORIGIN[which(ind.df$ORIGIN == "*")] = "I" # treat individuals with missing origin as immigrants (in the sense that they spontaneously appear in the population)

# check parent info if present
parents = FALSE
if ("SIRE_ID" %in% colnames(ind.df) == TRUE & "DAM_ID" %in% colnames(ind.df) == FALSE) stop("SIRE_ID detected without DAM_ID. Must give both parents or no parent Information at all.\n")
if ("DAM_ID" %in% colnames(ind.df) == TRUE & "SIRE_ID" %in% colnames(ind.df) == FALSE) stop("DAM_ID detected without SIRE_ID. Must give both parents or no parent information at all.\n")
if ("SIRE_ID" %in% colnames(ind.df) == TRUE) parents = TRUE

anc.idx = which(ind.df$ANCESTOR == 1)
if (length(anc.idx) > 1) {
	if (length(unique(ind.df$SEX[anc.idx])) < 2) stop("Error: Ancestral pairing impossible. Less than 2 sex provided.\n")
} else stop("Error: ind_file contains ", length(anc.idx), " ancestors. Poplation needs to be seeded with at least 2 ancestors\n")

# years to simulate
years = seq(from=ind.df$COHORT[1],to=ind.df$COHORT[1]+params$n_timesteps-1, by=1)

# dataframe of fledging probabilities per breeding pair

reproduction.df = fmtProbDist(read.table(params$offspring_file, head=FALSE, comment.char="#", sep="\t"))

# dataframe of lifespan probabilities

survive.df = NULL
if (!is.na(params$survive_file)) survive.df = fmtProbDist(read.table(params$survive_file, head=FALSE, comment.char="#", sep="\t"))

# parse events 

events = NULL
pair_limits = NA
maxpairs = NA # limits pairing after a certain number of pairs exists

if (!is.na(params$events_file)) {
	events = read.table(file=params$events_file, head=TRUE, comment.char="#", sep="\t")
	colnames(events) = toupper(colnames(events))
	if ('TIME' %in% colnames(events) == FALSE) stop("No TIME found in --events_file")
	idx = which(colnames(events) == 'TIME')
	idx_all = 1:ncol(events)
	events = events[,c(idx,idx_all[-which(idx_all == idx)])]
	events = events[order(events[,1]),]
	for (i in 1:ncol(events)) {
		if (colnames(events)[i] == 'TIME') {
			events[,i] = as.numeric(events[,i])
			if (length(which(events[,i] < 1)) > 0) stop("Error: --events_file contains TIME < 1\n")
		} else if (colnames(events)[i] == 'PAIR_LIMITS') {
			pair_limits = i
			events[,i] = toupper(events[,i])
			events[which(events[,i] == 'NA' || events[,i] == '*'),i] = NA
			events[,i] = as.numeric(events[,i])
		} else {
			stop(paste0("Unknown event ", colnames(events)[i], "in --events_file"))
		}
	}
}

## set up information storage containers

# ID maps

n_anc = length(which(ind.df$ANCESTOR == 1))
idmap = data.frame(matrix(NA, nrow=3*n_anc, ncol=2)) # maps IDs of individuals outside of or ancestral to the focal population to their ID in the simulation
idmap$X2 = as.numeric(idmap$X2)

if (parents) {
	ind.sub = ind.df[which((ind.df$ANCESTOR == 1 & ind.df$ORIGIN != "LR") | (ind.df$ANCESTOR == 0 & ind.df$ORIGIN == "I")),]
	external_parents = unique(c(ind.sub$SIRE_ID[which(ind.sub$SIRE_ID %in% ind.sub$ID == FALSE)], ind.sub$DAM_ID[which(ind.sub$DAM_ID %in% ind.sub$ID == FALSE)]))
	external_parents = external_parents[-which(external_parents == '*')]
	for (i in 1:length(external_parents)) {
		idmap$X1[i] = external_parents[i]
		idmap$X2[i] = i
	}
}

# survival distributions

survive.dist = list()
if (params$empiric_survive == TRUE) {
	for (i in 1:params$n_timesteps) {
		yr = ind.df$COHORT[1]+i-1
		idx = which(ind.df$COHORT == yr & ind.df$ORIGIN == "LR")
		if (length(idx) > 0) {
			survive.dist[[i]] = ind.df$COHORT_LAST[idx]
		} else {
			survive.dist[[i]] = NA
		}
	}
}
names(survive.dist) = years

## simulate population

pedinfo = c("ID", "SIRE_ID", "DAM_ID", "SEX", "COHORT", "COHORT_LAST", "YEAR_CLASS")
nr = 50
nc = length(pedinfo)
ped = data.frame(matrix(NA, nrow=length(which(ind.df$ANCESTOR == 1))+nr, ncol=nc)) # population pedigree
colnames(ped) = pedinfo
ped$ID = as.character(ped$ID)
ped$SIRE_ID = as.character(ped$SIRE_ID)
ped$DAM_ID = as.character(ped$DAM_ID)
ped$SEX = as.character(ped$SEX)
ped$COHORT = as.numeric(ped$COHORT)
ped$COHORT_LAST = as.numeric(ped$COHORT_LAST)
ped$YEAR_CLASS = as.numeric(ped$YEAR_CLASS)

lastmap = which(!is.na(idmap$X2))
simid = ifelse(length(lastmap) > 0, idmap$X2[lastmap[length(lastmap)]], 0)

breedpair = data.frame(matrix(NA, nrow=nr, ncol=3, dimnames=list(NULL,c("male","female","offspring")))) # matrix of breeding pairs [male ID, female ID, Number offspring]
breedpair$male = as.character(breedpair$male)
breedpair$female = as.character(breedpair$female)
breedpair$offspring = as.numeric(breedpair$offspring)
mature_pool = list(males = numeric(), females = numeric()) # set of mature available breeders

popstats = data.frame(year=years, n = NA, n_adult = NA, n_breeder = NA, n_pairs = NA) # note that n_breeders is the number of individuals actually contributing to next generation

for (time_step in 1:params$n_timesteps) {
	yr = years[time_step]
	
	#print(yr) # debug

	# insert non-local recruit ancestor, translocated, and immigrant individuals into the population

	idx = c(which(ind.df$COHORT == yr & ind.df$ANCESTOR == 1 & ind.df$ORIGIN != "LR"), which(ind.df$COHORT == yr & ind.df$ANCESTOR == 0 & ind.df$ORIGIN == "I"))
	if (length(idx) > 0) {
		# add IDs to the ID map if needed
		anc.id = ind.df$ID[which(ind.df$ID[idx] %in% ind.df$ID[which(ind.df$ANCESTOR == 1)] == TRUE)]
		idx2 = which(anc.id %in% idmap$X1 == FALSE)
		if (length(idx2) > 0) {
			start = which(is.na(idmap))[1]
			end = start + length(idx2) - 1
			idmap[start:end,1] = anc.id[idx2]
			idmap[start:end,2] = (simid+1):(simid+length(idx2))
			simid = simid+length(idx2)
		}

		# add individuals to pedigree
		if (length(idx) > length(which(is.na(ped$ID)))) ped = rbind(ped, matrix(NA, nrow=nr, ncol=nc, dimnames=list(NULL,colnames(ped)))) 
		ps = which(is.na(ped$ID))[1]
		pe = ps + length(idx) - 1
		yrlast = ind.df$COHORT_LAST[idx]
		yrlast[which(yrlast > years[length(years)])] = years[length(years)]
		ped[ps:pe,c(1,4,5,6,7)] = ind.df[idx,c("ID", "SEX", "COHORT", "COHORT_LAST", "YEAR_CLASS")]
		for (j in ps:pe) {
			mapidx = which(idmap$X1 == ped$ID[j])
			if (length(mapidx) > 0) {
				ped$ID[j] = idmap$X2[mapidx]
			} else {
				simid = simid+1
				ped$ID[j] = simid
			}
		}
		ped$SIRE_ID[ps:pe] = sapply(ind.df$SIRE_ID[idx], function(x,map){ifelse(x == '*', '*', map$X2[which(map$X1 == x)])}, map=idmap, USE.NAMES=FALSE)
		ped$DAM_ID[ps:pe] = sapply(ind.df$DAM_ID[idx], function(x,map){ifelse(x == '*', '*', map$X2[which(map$X1 == x)])}, map=idmap, USE.NAMES=FALSE)
	}
	ped$COHORT_LAST[which(ped$COHORT_LAST > years[length(years)])] = years[length(years)] # if individual lived longer than simulation, truncate at the last time point

	# update pool of available breeders
	mature_pool$males = ped$ID[which(ped$SEX == 'MALE' & (yr - ped$YEAR_CLASS) >= params$mature & ped$COHORT_LAST >= yr & ped$ID %in% breedpair[,1] == FALSE)]
	mature_pool$females = ped$ID[which(ped$SEX == 'FEMALE' & (yr - ped$YEAR_CLASS) >= params$mature & ped$COHORT_LAST >= yr & ped$ID %in% breedpair[,2] == FALSE)]
	
	#print(mature_pool) # debug

	# remove dead individuals
	rmv_id = ped$ID[which(ped$COHORT_LAST == (yr-1))] # these are individuals that are now dead in this timestep
	if (length(rmv_id) > 0) {
		# remove individuals from breeding pairs
		idx = which(breedpair[,1] %in% rmv_id == TRUE)
		if (length(idx) > 0) breedpair[idx,1] = NA
		idx = which(breedpair[,2] %in% rmv_id == TRUE)
		if (length(idx) > 0) breedpair[idx,2] = NA
		
		# return unpaired individuals to pool of available, unpaired breeders
		idx = which(!is.na(breedpair[,1]) & is.na(breedpair[,2]))
		if (length(idx) > 0) {
			mature_pool$males = c(mature_pool$males, as.vector(breedpair[idx,1]))
			breedpair[idx,1] = NA
		}

		idx = which(is.na(breedpair[,1]) & !is.na(breedpair[,2]))
		if (length(idx) > 0) {
			mature_pool$females = c(mature_pool$females, as.vector(breedpair[idx,2]))
			breedpair[idx,2] = NA
		}

		# remove individuals from pool of breeders
		idx = which(mature_pool$males %in% rmv_id == TRUE)
		if (length(idx) > 0) mature_pool$males = mature_pool$males[-idx]
		idx = which(mature_pool$females %in% rmv_id == TRUE)
		if (length(idx) > 0) mature_pool$females = mature_pool$females[-idx]
	}

	# form new breeding pairs
	nmature = length(mature_pool$males) + length(mature_pool$females)
	new_pairs = rpois(1,exp(params$pair_param[1]+params$pair_param[2]*nmature))
	if (!is.na(pair_limits) && time_step %in% events$TIME == TRUE) maxpairs = events[time_step,pair_limits] # set any potential new limits on forming breeding pairs
	breedpair[,3] = NA # for safety
	
	for (i in 1:new_pairs) {
		if (!is.na(maxpairs) && length(which(!is.na(breedpair[,1]))) >= maxpairs) break
		if (length(mature_pool$male) > 0 && length(mature_pool$female) > 0) {
			if (length(which(is.na(breedpair[,1]))) < new_pairs) breedpair = rbind(breedpair, matrix(NA, nrow=nr, ncol=3, dimnames=list(NULL,c("male","female","offspring"))))
			idx = which(is.na(breedpair[,1]))[1]
			breedpair[idx,1] = sample(mature_pool$males,size=1,replace=FALSE)
			mature_pool$males = mature_pool$males[-which(mature_pool$males == breedpair[idx,1])]
			breedpair[idx,2] = sample(mature_pool$females,size=1,replace=FALSE)
			mature_pool$females = mature_pool$females[-which(mature_pool$females == breedpair[idx,2])]
		} else break
	}

	# Reproduction
	
	# check if there are breeding pairs and determine the number of offspring for each breeding pair according to fledging number probabilities
	idx = which(!is.na(breedpair[,1])) # check if there are breeding pairs
	if (length(idx) > 0) {
		breedpair[idx,3] = sapply(runif(length(idx), min=0, max=1), function(x,df){if (x == 0) x = .Machine$double.eps; df$V1[which(x > df$a & x <= df$b)]}, df=reproduction.df, USE.NAMES=FALSE)
		overidx = which(breedpair[,3] > params$max_offspring)
		if (length(overidx) > 0) breedpair[overidx,3] = params$max_offspring
		noffspring = sum(breedpair[,3], na.rm=TRUE)

		# optionally redistribute offspring uniformly at random among breeding pairs
		if (params$unif_reproduction == TRUE) {
			maxidx = which(breedpair[,3] == params$max_offspring)
			if (length(maxidx) < length(idx)) {
				breedpair[idx,3] = 0
				assigned = 0
				while (assigned < noffspring) {
					availidx = which(breedpair[,3] < params$max_offspring)
					pairidx = as.numeric(sample(as.character(availidx),size=1,replace=FALSE))
					breedpair[pairidx,3] = breedpair[pairidx,3]+1
					assigned = assigned+1
				}
			}
		}
		
	}

	# update number of breeders contributing to the next generation info here (because the number of offspring is iterated later)
	popstats$n_breeder[time_step] = 2*length(which(!is.na(breedpair[,1]) & breedpair[,3] > 0))

	#print(breedpair) # debug

	# insert offspring into the pedigree

	# start by inserting local recruit ancestors
	lranc.idx = which(ind.df$COHORT == yr & ind.df$ANCESTOR == 1 & ind.df$ORIGIN == "LR") # determine if local recruits considered ancestors were born this timestep
	if (length(lranc.idx) > 0) {
		ancdf = ind.df[lranc.idx,]
		# try to assign siblings to the same nest if possible
		if (parents) ancdf$fam = paste0(ancdf$SIRE_ID,"_",ancdf$DAM_ID) else ancdf$fam = as.character(1:nrow(ancdf))
		famtab = sort(table(ancdf$fam),decreasing=TRUE)
		pairavail = which(!is.na(breedpair[,1])) # these are all pairs available to serve as parents
		
		if (length(pairavail) > 0) {
			# assign siblings  to pairs
			for (id in names(famtab)) {
				nsibs = unname(famtab[id]) 
				pairidx = which(breedpair[,3] >= nsibs) # all pairs with enough offspring to accommodate the siblings
				pairidx2 = NULL # index of breeding pair that will be the parents
				if (length(pairidx) > 0) {
					pairidx2 = as.numeric(sample(as.character(pairidx),size=1,replace=FALSE))
				} else {
					# need to artificially inflate the largest family size to accommodate all siblings
					pairidx2 = as.numeric(sample(as.character(which(breedpair[,3] == max(breedpair[,3],na.rm=TRUE))), size=1, replace=FALSE))
					breedpair[pairidx2,3] = nsibs
					if (breedpair[pairidx2,3] > params$max_offspring) {
						write("Warning: Number offspring increased beyond max_offspring to accommodate ancestral LR siblings",stderr())
					}
				}
				# insert LR ancecstors with known parents into pedigree
				for (k in which(ancdf$fam == id)) {
					if (length(which(is.na(ped$ID))) < 1) ped = rbind(ped, matrix(NA, nrow=nr, ncol=nc, dimnames=list(NULL,colnames(ped))))
					pedidx = which(is.na(ped$ID))[1]
					simid = simid+1
					ped$ID[pedidx] = simid
					ped[pedidx,2:3] = breedpair[pairidx2,1:2]
					ped[pedidx,4:7] = ancdf[k,c("SEX", "COHORT", "COHORT_LAST", "YEAR_CLASS")]
					idmap[which(is.na(idmap$X1))[1],] = c(ancdf$ID[k], simid)
					breedpair[pairidx2,3] = breedpair[pairidx2,3]-1
				}
				ancdf = ancdf[-which(ancdf$fam == id),]
				pairavail = pairavail[-which(pairavail == pairidx2)]
				if (length(pairavail) == 0) break # no breeders are available
			}
		}

		if (nrow(ancdf) > 0) {
			# LR ancestors insert but no breeders are available
			# have to migrate the LR ancestors into the population with missing parent information
			for (k in 1:nrow(ancdf)) {
				if (length(which(is.na(ped$ID))) < 1) ped = rbind(ped, matrix(NA, nrow=nr, ncol=nc, dimnames=list(NULL,colnames(ped))))
				pedidx = which(is.na(ped$ID))[1]
				simid = simid+1
				ped$ID[pedidx] = simid
				ped[pedidx,2:3] = c('*', '*')
				ped[pedidx,4:7] = ancdf[k,c("SEX", "COHORT", "COHORT_LAST", "YEAR_CLASS")]
				idmap[which(is.na(idmap$X1))[1],] = c(ancdf$ID[k], simid)
			}
		}

	}

	# insert remaining offspring into pedigree
	idx = which(breedpair[,3] > 0) # breeding pairs with offspring (that haven't yet been inserted into the pedigree)
	if (length(idx) > 0) {
		for (pairidx in idx) {
			if (length(which(is.na(ped$ID))) < breedpair[pairidx,3]) ped = rbind(ped, matrix(NA, nrow=nr, ncol=nc, dimnames=list(NULL,colnames(ped))))
			start = which(is.na(ped$ID))[1]
			end = start + breedpair[pairidx,3] - 1
			ped$ID[start:end] = (simid+1):(simid + breedpair[pairidx,3])
			ped$SIRE_ID[start:end] = breedpair[pairidx,1]
			ped$DAM_ID[start:end] = breedpair[pairidx,2]
			ped$SEX[start:end] = sample(c('MALE','FEMALE'),size=breedpair[pairidx,3],replace=TRUE, prob=c(params$male_p, 1-params$male_p))
			ped$COHORT[start:end] = yr
			if (params$empiric_survive == TRUE) {
				yridx = which(names(survive.dist) == yr)
				if (length(which(!is.na(survive.dist[[yridx]]))) >= params$min_survive_obs) {
					ped$COHORT_LAST[start:end] = sample(as.character(survive.dist[[yridx]]), size=breedpair[pairidx,3], replace=TRUE)
				} else {
					if (is.null(survive.df)) stop(cat("No survival information for ",yr," local recruits\n"))
					ped$COHORT_LAST[start:end] = yr + sapply(runif(breedpair[pairidx,3], min=0, max=1), 
					function(x,df){if (x == 0) x = .Machine$double.eps; df$V1[which(x > df$a & x <= df$b)]}, df=survive.df, USE.NAMES=FALSE)
				}
			} else {
				ped$COHORT_LAST[start:end] = yr + sapply(runif(breedpair[pairidx,3], min=0, max=1), 
				function(x,df){if (x == 0) x = .Machine$double.eps; df$V1[which(x > df$a & x <= df$b)]}, df=survive.df, USE.NAMES=FALSE)
			}
			ped$YEAR_CLASS[start:end] = yr
			simid = simid + breedpair[pairidx,3]
		}
	}

	#print(idmap) # debug
	#print(ped) # debug

	# update population statistics
	popstats$n[time_step] = length(which(ped$COHORT_LAST >= yr)) # total number of alive individuals
	popstats$n_pairs[time_step] = length(which(!is.na(breedpair[,1]))) # number of breeding pairs
	popstats$n_adult[time_step] = length(mature_pool$males) + length(mature_pool$females) + 2*popstats$n_pairs[time_step] # number of adults (individuals with age >= --mature)
}

# if individual lived longer than simulation, truncate last observation at the last time point
truncate_idx = which(ped$COHORT_LAST > years[length(years)])
if (length(truncate_idx) > 0) ped$COHORT_LAST[which(ped$COHORT_LAST > years[length(years)])] = years[length(years)]

write("\nSimulation complete", stderr())

## write output

write(paste0("\nWriting result files:\n", out_ped, "\n", out_stats, "\n", out_idmap), stderr())

write.table(ped[which(!is.na(ped$COHORT)),c(1:6)], file=out_ped, col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")
write.table(idmap[which(idmap$X1 %in% ind.df$ID[which(ind.df$ANCESTOR == 1)]),c(2,1)], file=out_idmap, col.names=FALSE, row.names=FALSE, quote=FALSE, sep="\t")
write.table(popstats, file=out_stats, col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")
