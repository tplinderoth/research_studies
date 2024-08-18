#!/bin/bash
#! This line is a comment
#! Make sure you only have comments and #SBATCH directives between here and the end of the #SBATCH directives, or things will break
#! Name of the job:
#SBATCH -J fsj_H_froh_ps_baw
#! Account name for group, use SL2 for paying queue:
#SBATCH -A general
#! Output filename:
#! %A means slurm job ID and %a means array index
#SBATCH --output=fsj_H_froh_ps_baq_%A_%a.out
#! Errors filename:
#SBATCH --error=fsj_H_froh_ps_baq_%A_%a.err

#! Number of nodes to be allocated for the job (for single core jobs always leave this at 1)
#SBATCH --nodes=1
#! Number of tasks. By default SLURM assumes 1 task per node and 1 CPU per task.
#SBATCH --ntasks=1
#! How many many cores will be allocated per task?
#SBATCH --cpus-per-task=12
#! Estimated runtime: hh:mm:ss (job is force-stopped after if exceeded):
#SBATCH --time=06:00:00
#! Estimated maximum memory needed (job is force-stopped if exceeded):
#SBATCH --mem=48000mb
#! Submit a job array with index values between 0 and 31
#! NOTE: This must be a range, not a single number
#SBATCH --array=1-87

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
#! module load Java/1.8.0_152 bzip2/1.0.6 zlib/1.2.11 Boost/1.67.0 GSL/2.6

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

BAMLIST='/mnt/research/Fitz_Lab/projects/mosaic/map/mosaic_bam_list.txt'
BAMFILE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$BAMLIST")
SAMPLIST='/mnt/research/Fitz_Lab/projects/mosaic/mosiac_ids.txt'
SAMP=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$SAMPLIST")
REF='/mnt/research/Fitz_Lab/ref/bird/FSJ_V3/FSJ.V3.fa'
REGFILE='/mnt/research/Fitz_Lab/ref/bird/FSJ_V3/FSJ_V3_main_autosomes.txt'
SITESFILE='/mnt/research/Fitz_Lab/projects/mosaic/variants/vcf/all_sites/fsj_mosaic_allsites_main_autosomes_qc.pos'
OUTPREFIX="/mnt/research/Fitz_Lab/projects/mosaic/popgen/theta/heterozygosity/Froh_H/sfs/pop_specific/${SAMP}_allsites_main_autosomes_qc_Froh_pop_specific_baq"
FFILE="/mnt/research/Fitz_Lab/projects/mosaic/popgen/fstats/individual_f_files/${SAMP}_biallelic_snps_main_autosomes_qc_pop_specific_af_Froh.txt"
EXEC='/mnt/research/Fitz_Lab/projects/mosaic/popgen/theta/heterozygosity/ngsF_H/scripts/saf2H.pl'
EXEC2='/mnt/research/Fitz_Lab/projects/mosaic/popgen/sfs/revision1_sfs/saf2sfs.pl'

SAFCMD="angsd -i $BAMFILE -baq 1 -C 50 -ref $REF -out $OUTPREFIX -GL 1 -doSaf 2 -indF $FFILE -anc $REF -doMaf 1 -doMajorMinor 1 -minQ 20 -minMapQ 20 -remove_bads 1 -only_proper_pairs 1 -rf $REGFILE -sites $SITESFILE -P 4"

printf "\n%s\n\n" "$SAFCMD"

eval $SAFCMD
wait

HCMD="realSFS print ${OUTPREFIX}.saf.idx | $EXEC > ${OUTPREFIX}_H.txt"

printf "\n%s\n\n" "$HCMD"

eval $HCMD
wait

SFSCMD="realSFS print ${OUTPREFIX}.saf.idx | $EXEC2 0 > ${OUTPREFIX}_unfold.sfs"

printf "\n%s\n\n" "$SFSCMD"

eval $SFSCMD

