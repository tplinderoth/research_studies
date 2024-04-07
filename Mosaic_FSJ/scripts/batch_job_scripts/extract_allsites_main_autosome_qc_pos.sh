#!/bin/bash
#! This line is a comment
#! Make sure you only have comments and #SBATCH directives between here and the end of the #SBATCH directives, or things will break
#! Name of the job:
#SBATCH -J allsites_autosomes_qc_pos
#! Account name for group, use SL2 for paying queue:
#! #SBATCH -A ACCOUNT_NAME

#! Number of nodes to be allocated for the job (for single core jobs always leave this at 1)
#SBATCH --nodes=1
#! Number of tasks. By default SLURM assumes 1 task per node and 1 CPU per task. (for single core jobs always leave this at 1)
#SBATCH --ntasks=1
#! How many many cores will be allocated per task? (for single core jobs always leave this at 1)
#SBATCH --cpus-per-task=6
#! Estimated runtime: hh:mm:ss (job is force-stopped after if exceeded):
#SBATCH --time=36:00:00
#! Estimated maximum memory needed (job is force-stopped if exceeded):
#! #SBATCH --mem=66000mb
#! #SBATCH --no-requeue

#! This is the partition name.
#! #SBATCH -p skylake-himem

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
#! safe value to no more than 32:
#! export OMP_NUM_THREADS=1

echo "This is job ${SLURM_JOB_ID}"

workdir="$SLURM_SUBMIT_DIR" # The value of SLURM_SUBMIT_DIR sets workdir to the directory
cd $workdir

VCF='/mnt/research/Fitz_Lab/projects/mosaic/variants/vcf/all_sites/fsj_mosaic_allsites_genome.vcf.gz'
OUTPOS='/mnt/research/Fitz_Lab/projects/mosaic/variants/vcf/all_sites/fsj_mosaic_allsites_main_autosomes_qc.pos'
REGFILE='/mnt/research/Fitz_Lab/ref/bird/FSJ_V3/fsj_v3_main_autosomes.rf'

bcftools view -H -R "$REGFILE" -f "PASS" -i 'N_PASS(GT[0-27]!="mis" & FMT/DP > 2) > 24 && N_PASS(GT[28-38]!="mis" & FMT/DP > 2) > 9 && N_PASS(GT[39-44]!="mis" & FMT/DP > 2) > 4 && N_PASS(GT[45-57]!="mis" & FMT/DP > 2) > 11 && N_PASS(GT[58-86]!="mis" & FMT/DP > 2) > 25' "$VCF" | perl -ne  '@tok = split(/\s+/,$_); if ($tok[7] =~ /REPGQ=(\d+),(\d+),(\d+),(\d+),(\d+)/) {print "$tok[0]\t$tok[1]\n" if ($1 > 24 && $2 > 9 && $3 > 4 && $4 > 11 && $5 > 25);} else {print "$tok[0]\t$tok[1]\n";}' | uniq > "$OUTPOS"

# Site extraction command breakdown:

# take only sites with FILTER=PASS and for which 90% of individuals in each C,E,I,M4,T class have non-missing genotypes and minimum DP of 3
# bcftools view -H -R "$REGFILE" -f "PASS" -i 'N_PASS(GT[0-27]!="mis" & FMT/DP > 2) > 24 && N_PASS(GT[28-38]!="mis" & FMT/DP > 2) > 9 && N_PASS(GT[39-44]!="mis" & FMT/DP > 2) > 4 && N_PASS(GT[45-57]!="mis" & FMT/DP > 2) > 11 && N_PASS(GT[58-86]!="mis" & FMT/DP > 2) > 25' "$VCF"

# If INFO/REPGQ flag is present require 90% of individuals in each C,E,I,M4,T to have a minimum genotype quality of 15 in order to retain site
# perl -ne  '@tok = split(/\s+/,$_); if ($tok[7] =~ /REPGQ=(\d+),(\d+),(\d+),(\d+),(\d+)/) {print "$tok[0]\t$tok[1]\n" if ($1 > 24 && $2 > 9 && $3 > 4 && $4 > 11 && $5 > 25);} else {print "$tok[0]\t$tok[1]\n";}

# uniq is because some sites show up more than once in the VCF

angsd sites index "$OUTPOS"
