#!/bin/bash
#! This line is a comment
#! Make sure you only have comments and #SBATCH directives between here and the end of the #SBATCH directives, or things will break
#! Name of the job:
#SBATCH -J sim_pedigree_cohort_contributions
#! Account name for group, use SL2 for paying queue:
#! #SBATCH -A ACCOUNT_NAME
#! Output filename:
#! %A means slurm job ID and %a means array index
#SBATCH --output=sim_ped_cohort_analysis_%A_%a.out
#! Errors filename:
#SBATCH --error=sim_ped_cohort_analysis_%A_%a.err

#! Number of nodes to be allocated for the job (for single core jobs always leave this at 1)
#SBATCH --nodes=1
#! Number of tasks. By default SLURM assumes 1 task per node and 1 CPU per task.
#SBATCH --ntasks=1
#! How many many cores will be allocated per task?
#SBATCH --cpus-per-task=8
#! Estimated runtime: hh:mm:ss (job is force-stopped after if exceeded):
#SBATCH --time=36:00:00
#! Estimated maximum memory needed (job is force-stopped if exceeded):
#SBATCH --mem=16000mb
#! Submit a job array with index values between 0 and 31
#! NOTE: This must be a range, not a single number
#SBATCH --array=1-10

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

SEED_F="/mnt/research/Fitz_Lab/projects/mosaic/simulations/pedigree_sims/inputs/seeds/seeds_${SLURM_ARRAY_TASK_ID}.txt"
PED_F='/mnt/research/Fitz_Lab/projects/mosaic/simulations/pedigree_sims/inputs/cr_ped_20240602.tsv'
JAY_META_F='/mnt/research/Fitz_Lab/projects/mosaic/simulations/pedigree_sims/inputs/jayid_metadata_20240602.tsv'
NEST_META_F='/mnt/research/Fitz_Lab/projects/mosaic/simulations/pedigree_sims/inputs/nest_metadata_20240602.tsv'
ANC_F='/mnt/research/Fitz_Lab/projects/mosaic/simulations/pedigree_sims/inputs/cr_ancestors_20240602.tsv'
FOCAL_F='/mnt/research/Fitz_Lab/projects/mosaic/simulations/pedigree_sims/inputs/res_trans_focal_IDs.txt'
SIMPED='/mnt/research/Fitz_Lab/projects/mosaic/simulations/pedigree_sims/scripts/simPed.R'
RELSTATS='/mnt/research/Fitz_Lab/software/PopGenomicsTools/relateStats'
SIM_WRAPPER='/mnt/research/Fitz_Lab/projects/mosaic/simulations/pedigree_sims/scripts/pedSimStats_v1.pl'
OUTPREFIX="/mnt/research/Fitz_Lab/projects/mosaic/simulations/pedigree_sims/results/cr_sim_cohort_analysis_${SLURM_ARRAY_TASK_ID}"

CMD="$SIM_WRAPPER \
--ped_file $PED_F \
--indmeta $JAY_META_F \
--nestmeta $NEST_META_F \
--ancfile $ANC_F \
--out $OUTPREFIX \
--simexec $SIMPED \
--statsexec $RELSTATS \
--focal_ind $FOCAL_F \
--seed_file $SEED_F \
--rp 0 \
--maxoffspring 5 \
--mature 2 \
--p_male 0.5 \
--cohort \
2>/dev/null"

printf "\n%s\n\n" "$CMD"

eval $CMD
