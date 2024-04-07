#!/bin/bash

# map_reads.sh <reference fasta> <sample name, e.g. 'T_101'> <fastq directory> <output directory> [number threads]

if [ $# -lt 4 ]
then
	>&2 echo "map_reads.sh <reference fasta> <sample name, e.g. 'T_101'> <fastq directory> <output directory> [number threads]"
	exit
fi

NTHREAD=1
REF=$1 # reference fasta
SAMPLE_NAME=$2 # sample name, e.g. T_101
FQDIR=$3
OUTDIR=$4
if [ $# -gt 4 ]; then NTHREAD=$5; fi

DATADS='Whole genome sequences from Florida scrub jays sequenced for the Mosaic project.'

if [[ "$FQDIR" == */ ]]; then FQDIR=$(echo "$FQDIR" | sed 's/\/$//'); fi

if [[ "$OUTDIR" == */ ]]; then OUTDIR=$(echo "$OUTDIR" | sed 's/\/$//'); fi
OUTPREFIX="${OUTDIR}/${SAMPLE_NAME}"

# get read group info from fastq files

fqlist=($(find $FQDIR -name "${SAMPLE_NAME}*.fastq.gz" | tr '\n' ' '))
numlanes=$(perl -e '$nfile = <>; print $nfile/2' <<< ${#fqlist[@]}) # this is the nubmer of lanes library was sequenced on
if [ $numlanes -lt 1 ]; then >&2 printf "No files with sample name %s in %s\n" "$SAMPLE_NAME" "$FQDIR"; exit 1; fi
samplib=$(echo "${fqlist[0]}" | perl -e '$fqname = <>; print $1 if $fqname =~ /_([^_]+)_\d+_[^_]+\.fastq\.gz$/')
fqprefix="${SAMPLE_NAME}_${samplib}"
declare -a readheader

# LOOP OVER LANES TO MAP EACH SET OF PE FILES
for lanenum in $(seq 1 $numlanes)
do

	prevrg='N'
	readgroup=''
	nmatch=0
	lstart=-3

	# need to check multiple headers to make sure sample barcodes converage (i.e. no sequencing error)
	while [ $nmatch -lt 10 ]
	do
		lstart=$((lstart + 4))
		lend=$((lstart + 1))
		readheader=$(zcat "${FQDIR}/${fqprefix}_${lanenum}_R1.fastq.gz" | sed -n "${lstart}p;${lend}q" | sed 's/ [[:digit:]]//')
		readinfo=($(echo "$readheader" | perl -e 'chomp($read = <>); @arr = split(/:/,$read); print "@arr[2,3,$#arr]";'))
		barcode=$(sed 's/+/-/' <<< "${readinfo[2]}")
		readgroup="@RG\tID:${readinfo[0]}.${readinfo[1]}\tBC:${barcode}\tDS:${DATADS}\tLB:${samplib}\tPL:ILLUMINA\tPU:${readinfo[0]}.${readinfo[1]}.${barcode}\tSM:${SAMPLE_NAME}"

		if [ "$readgroup" = "$prevrg" ] && ! [[ "$barcode" =~ 'N' ]]; then ((nmatch++)); fi
		prevrg="$readgroup"
	done

	printf "\n--MAPPING READS (lane %i/%i)--\n" "$lanenum" "$numlanes" 

	# map paired-end reads

	FWDFQ="${FQDIR}/${SAMPLE_NAME}_${samplib}_${lanenum}_R1.fastq.gz"
	REVFQ="${FQDIR}/${SAMPLE_NAME}_${samplib}_${lanenum}_R2.fastq.gz"

	OUTBAM="${OUTDIR}/${SAMPLE_NAME}_raw"
	if [ $numlanes -gt 1 ]; then OUTBAM+="_${lanenum}.bam"; else OUTBAM+=".bam"; fi

	MAPCMD="bwa mem -t $NTHREAD -R '$readgroup' $REF $FWDFQ $REVFQ | samtools view -b > $OUTBAM"

	printf "\n%s\n\n" "$MAPCMD"

	eval $MAPCMD

	wait

done

printf "\n--MERGING SEPARATE LANE BAMS--\n"

if [ $numlanes -gt 1 ]; then
	declare -a mergelist
	for i in $(seq 1 $numlanes)
	do
		mergelist+=("${OUTDIR}/${SAMPLE_NAME}_raw_${i}.bam")
	done
	MERGE_CMD="samtools merge -n -c -@ $NTHREAD -o ${OUTDIR}/${SAMPLE_NAME}_raw.bam ${mergelist[@]}"
	printf "\n%s\n\n" "$MERGE_CMD"
	eval $MERGE_CMD
else
	printf "\nSingle lane, no merging necessary.\n"
fi

printf "\n--SORTING AND MARKING DUPLICATES--\n"

# The following pipe
# 1) Adds ms (mate score) tag for duplicate marking (bwa output should already be grouped by read name)
# 2) sorts the bam
# 3) marks duplicates

DIVTHREAD=$(perl -e '$thread = <>; print int($thread/3)' <<< "$NTHREAD")
OUTBAM2="${OUTDIR}/${SAMPLE_NAME}_untrimmed.bam"

POSTMAP_CMD="samtools fixmate -m -@ $DIVTHREAD -u ${OUTDIR}/${SAMPLE_NAME}_raw.bam - | samtools sort -u -@ $DIVTHREAD | samtools markdup -s -S -O BAM -@ $DIVTHREAD - ${OUTBAM2}"

printf "\n%s\n\n" "$POSTMAP_CMD"

eval $POSTMAP_CMD

printf "\n--INDEXING BAM--\n"

INDEX_CMD="samtools index -@ $NTHREAD ${OUTBAM2}"

printf "\n%s\n" "$INDEX_CMD"

eval $INDEX_CMD

printf "\n--FINISHED--\n"

exit
