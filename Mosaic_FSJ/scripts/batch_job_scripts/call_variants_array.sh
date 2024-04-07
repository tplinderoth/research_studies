#!/bin/bash
#! This line is a comment
#! Make sure you only have comments and #SBATCH directives between here and the end of the #SBATCH directives, or things will break
#! Name of the job:
#SBATCH -J mosaic_variant_call
#! Account name for group, use SL2 for paying queue:
#! #SBATCH -A ACCOUNT_NAME
#! Output filename:
#! %A means slurm job ID and %a means array index
#SBATCH --output=mosaic_variant_call_%A_%a.out
#! Errors filename:
#SBATCH --error=mosaic_variant_call_%A_%a.err

#! Number of nodes to be allocated for the job (for single core jobs always leave this at 1)
#SBATCH --nodes=1
#! Number of tasks. By default SLURM assumes 1 task per node and 1 CPU per task.
#SBATCH --ntasks=1
#! How many many cores will be allocated per task?
#SBATCH --cpus-per-task=24
#! Estimated runtime: hh:mm:ss (job is force-stopped after if exceeded):
#SBATCH --time=48:00:00
#! Estimated maximum memory needed (job is force-stopped if exceeded):
#SBATCH --mem=66000mb
#! Submit a job array with index values between 0 and 31
#! NOTE: This must be a range, not a single number
#SBATCH --array=1-15

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

REF='/mnt/research/Fitz_Lab/ref/bird/FSJ_V3/FSJ.V3.fa'
BAMS='/mnt/research/Fitz_Lab/projects/mosaic/map/mosaic_bam_list.txt'
SCAFFLIST="/mnt/research/Fitz_Lab/projects/mosaic/variants/call/scaffold_sets/regions_${SLURM_ARRAY_TASK_ID}.rf"
PLOIDY_FILE='/mnt/research/Fitz_Lab/projects/mosaic/variants/call/mosaic_ploidy.txt'
SAMP_FILE='/mnt/research/Fitz_Lab/projects/mosaic/variants/call/mosaic_sample_sex.txt'
OUTBCF="/mnt/gs18/scratch/users/lindero1/mosaic/vcf/call/fsj_mosaic_${SLURM_ARRAY_TASK_ID}.bcf.gz"

CMD="bcftools mpileup \
-f $REF \
-b $BAMS \
-R $SCAFFLIST \
-C 0 \
-d 10000 \
-L 10000 \
-q 20 \
-Q 13 \
--ns UNMAP,SECONDARY,QCFAIL,DUP \
-a FORMAT/AD,FORMAT/DP,FORMAT/QS,FORMAT/SP,FORMAT/SCR,INFO/AD,INFO/SCR \
-p \
-O u \
| bcftools call \
--ploidy-file $PLOIDY_FILE \
-S $SAMP_FILE \
-a PV4,GQ,GP \
-m \
-P 0.003 \
-O u \
| bcftools +fill-tags \
-O b \
-o $OUTBCF \
-- -t 'AF,ExcHet,NS'"

printf "\n%s\n\n" "$CMD"

eval $CMD

wait

tabix -p bcf $OUTBCF
