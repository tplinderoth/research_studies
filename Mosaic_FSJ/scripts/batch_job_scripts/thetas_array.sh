#!/bin/bash
#! This line is a comment
#! Make sure you only have comments and #SBATCH directives between here and the end of the #SBATCH directives, or things will break
#! Name of the job:
#SBATCH -J mosaic_thetas
#! Account name for group, use SL2 for paying queue:
#! #SBATCH -A ACCOUNT_NAME
#! Output filename:
#! %A means slurm job ID and %a means array index
#SBATCH --output=mosaic_thetas_%A_%a.out
#! Errors filename:
#SBATCH --error=mosaic_thetas_%A_%a.err

#! Number of nodes to be allocated for the job (for single core jobs always leave this at 1)
#SBATCH --nodes=1
#! Number of tasks. By default SLURM assumes 1 task per node and 1 CPU per task.
#SBATCH --ntasks=1
#! How many many cores will be allocated per task?
#SBATCH --cpus-per-task=24
#! Estimated runtime: hh:mm:ss (job is force-stopped after if exceeded):
#SBATCH --time=18:00:00
#! Estimated maximum memory needed (job is force-stopped if exceeded):
#SBATCH --mem=72000mb
#! Submit a job array with index values between 0 and 31
#! NOTE: This must be a range, not a single number
#SBATCH --array=0-6

#! This is the partition name.
#! #SBATCH -p cclake

#! mail alert at start, end and abortion of execution
#! emails will default to going to your email address
#! you can specify a different email address manually if needed.
#SBATCH --mail-type=FAIL

#! Don't put any #SBATCH directives below this line

#! Modify the environment seen by the application. For this example we need the default modules.
#! . /etc/profile.d/modules.sh                # This line enables the module command
#! module purge                               # Removes all modules still loaded
#! module load rhel7/default-peta4            # REQUIRED - loads the basic environment
module load Java/1.8.0_152 bzip2/1.0.6 zlib/1.2.11 Boost/1.67.0 GSL/2.6

#! Are you using OpenMP (NB this is unrelated to OpenMPI)? If so increase this
#! export OMP_NUM_THREADS=1

#! The variable $SLURM_ARRAY_TASK_ID contains the array index for each job.
#! In this example, each job will be passed its index, so each output file will contain a different value
echo "This is job" $SLURM_ARRAY_TASK_ID

#! Command line that we want to run:
#! jobDir=Job_$SLURM_ARRAY_TASK_ID
#! mkdir $jobDir
#! cd $jobDir

workdir="$SLURM_SUBMIT_DIR" # The value of SLURM_SUBMIT_DIR sets workdir to the directory
cd $workdir

GRP_ARR=('C' 'contemporary_unrelated' 'resident' 'resident_unrelated' 'T' 'translocated_unrelated' 'M4')
FSJGRP="${GRP_ARR[$SLURM_ARRAY_TASK_ID]}"
SAFIDX="/mnt/research/Fitz_Lab/projects/mosaic/popgen/sfs/fsj_mosaic_allsites_main_autosomes_qc_${FSJGRP}.saf.idx"
SFS="/mnt/research/Fitz_Lab/projects/mosaic/popgen/sfs/fsj_mosaic_allsites_main_autosomes_qc_${FSJGRP}.fold.sfs"
OUTFILE="/mnt/research/Fitz_Lab/projects/mosaic/popgen/theta/fsj_mosaic_allsites_main_autosomes_qc_${FSJGRP}"
THETAIDX="${OUTFILE}.thetas.idx"

CMD="realSFS saf2theta $SAFIDX -outname $OUTFILE -sfs $SFS -fold 1"

printf "\n%s\n\n" "$CMD"

eval $CMD
wait

CMD2="thetaStat do_stat $THETAIDX"

printf "\n%s\n\n" "$CMD2"

eval $CMD2
