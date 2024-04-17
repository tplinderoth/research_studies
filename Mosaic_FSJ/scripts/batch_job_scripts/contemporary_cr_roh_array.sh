#!/bin/bash
#! This line is a comment
#! Make sure you only have comments and #SBATCH directives between here and the end of the #SBATCH directives, or things will break
#! Name of the job:
#SBATCH -J contemporary_cr_roh_contemporary_cr_af
#! Account name for group, use SL2 for paying queue:
#! #SBATCH -A ACCOUNT_NAME
#! Output filename:
#! %A means slurm job ID and %a means array index
#SBATCH --output=contemporary_cr_roh_contemporary_cr_af_%A_%a.out
#! Errors filename:
#SBATCH --error=contemporary_cr_roh_contemporary_cr_af_%A_%a.err

#! Number of nodes to be allocated for the job (for single core jobs always leave this at 1)
#SBATCH --nodes=1
#! Number of tasks. By default SLURM assumes 1 task per node and 1 CPU per task.
#SBATCH --ntasks=1
#! How many many cores will be allocated per task?
#SBATCH --cpus-per-task=16
#! Estimated runtime: hh:mm:ss (job is force-stopped after if exceeded):
#SBATCH --time=36:00:00
#! Estimated maximum memory needed (job is force-stopped if exceeded):
#SBATCH --mem=72000mb
#! Submit a job array with index values between 0 and 31
#! NOTE: This must be a range, not a single number
#SBATCH --array=1-28

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

# recombination rate estimates for birds:
# Backstrom etal 2010: mean 1.5 cM/Mb in zebra finch exluding 10 microchromsomes and projected to be ~2 cM/Mb with the microchromosomes
# Singhal etal 2015: median 0.14 cM/Mb in both Taeniopygia guttata (zebra finch) and Poephila acuticuada (long-tailed finch)
# Bascon-Cardozo etal 2024: mean 5.8 cM/Mb in Eurasian blackcap
# Smeds etal 2016: mean 3.08 cM/Mb in collared flycatcher
# Ellegren 2005: median 2.8 cM/Mb for macrochromosomes in chicken
# Hagen etal 2020: mean 1.78 cM/Mb for macrochromosomes and 6.41 cM/Mb for microchromosomes, genome average = 4.095 cM/Mb, in house sparrows
# Groenen etal 2000: mean 6.02 cM/Mb in chicken
# Groenen etal 2009: mean 3.11 cM/Mb in chicken
# Stapley etal 2008: mean 3.18 cM/Mb in zebra finch
# Kawakami etal 2014: mean 3.1 cM/Mb in collared flycatcher

# Used in average estimate
# Backstrom etal 2010: 2 cM/Mb in zebrafinch
# Bascon-Cardozo etal 2024: 5.8 cM/Mb in Eurasian blackcap
# Smeds etal 2016: 3.08 cM/Mb in collared flycatcher
# Hagen etal 2020: 4.095 in house sparrows
# Groenen etal 2009: 3.11 cM/Mb in chicken

# Will use average of the mean recombination rates for the different species above:
# mean(c(2, 5.8, 3.08, 4.095, 3.11)) == 3.617 cM/Mb == 3.6e-8 crossovers/bp (expected)

VCF='/mnt/research/Fitz_Lab/projects/mosaic/variants/vcf/biallelic_snps/fsj_mosaic_biallelic_snps_main_autosomes_qc.vcf.gz'
SAMPLE_LIST='/mnt/research/Fitz_Lab/projects/mosaic/popgen/roh/inputs/contemporary_cr_ids.txt'
SAMP=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$SAMPLE_LIST")
AF_FILE='/mnt/research/Fitz_Lab/projects/mosaic/popgen/roh/inputs/contemporary_cr_maxr0.4_biallelic_snps_main_autosomes_qc.af.gz'
OUTFILE="/mnt/research/Fitz_Lab/projects/mosaic/popgen/roh/estimates/${SAMP}_contemporary_cr_af_roh.txt"

CMD="bcftools roh --samples $SAMP --AF-file $AF_FILE --rec-rate 3.6e-8 --output $OUTFILE --output-type sr --skip-indels --viterbi-training 1e-10 --threads 12 $VCF"

printf "\n%s\n\n" "$CMD"

eval $CMD
