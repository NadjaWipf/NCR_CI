#!/usr/bin/bash

SAMPLE=`sed -n ${SLURM_ARRAY_TASK_ID}p  sampleNames.txt`

InFOLDER="bamfiles/bam_00"
OutFOLDER="variant_01_merged_bam_files"

mkdir -p $OutFOLDER

ml SAMtools

samtools merge -X $OutFOLDER/$SAMPLE.merged.Aligned.sortedByCoord.out.bam ${InFOLDER}1/${SAMPLE}_001Aligned.sortedByCoord.out.bam ${InFOLDER}2/${SAMPLE}_002Aligned.sortedByCoord.out.bam

