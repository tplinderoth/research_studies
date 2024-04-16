#!/usr/bin/env Rscript

# parse inputs
args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 2) {
	cat("\nsubsetUnrelated.R <relatedness matrix with sample ID column names> <max relatedness> [seed]\n\n")
	options(show.error.messages=FALSE)
	stop()
}

rmat <- as.matrix(read.table(args[1], head=TRUE))
rmax <- as.numeric(args[2])
seed = sample(1:5000,1,replace=FALSE)
if (length(args) > 2) {
	seed = as.numeric(args[3])
	set.seed(seed)
}
write(paste0("random seed: ",seed),stderr())

# mask diagonal
diag(rmat) <- NA

submat.keep = NULL
for (k in 1:100) {
	submat = rmat
	while (nrow(submat) > 1 && max(submat,na.rm=TRUE) > rmax) {
		reset = 0
		#for (i in 1:nrow(submat)) {
		for (i in sample(1:nrow(submat),replace=FALSE)) {
			#for (j in 1:ncol(submat)) {
			for (j in sample(1:ncol(submat),replace=FALSE)) {
				if (i == j) next
				if (submat[i,j] > rmax) {
					rmv.idx = ifelse(length(which(submat[i,] > rmax)) > length(which(submat[j,] > rmax)), i, j)
					submat = submat[-rmv.idx, -rmv.idx]
					reset = 1;
					break;
				}
			}
			if (reset == 1) break;
		}
	}
	if (is.null(submat.keep) || nrow(submat) > nrow(submat.keep)) submat.keep = submat
}

write(colnames(submat.keep), stdout())
