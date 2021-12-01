#!/usr/bin/bash

ml SAMtools
ml picard

SAMPLE=`sed -n ${SLURM_ARRAY_TASK_ID}p  sampleNames.txt`

InFOLDER="variant_01_merged_bam_files"
OutFOLDER="variant_02_bamPreProcessing"

InBAM=$SAMPLE.merged.Aligned.sortedByCoord.out.bam

#############################################################################################
# PICARD - add read groups, sorting, remove duplicates

# Picard AddOrReplaceReadGroups & sorting
java -jar $EBROOTPICARD/picard.jar AddOrReplaceReadGroups I=$InFOLDER/$InBAM O=$OutFOLDER/$InBAM.rg.sort.bam \
		SO=coordinate RGID=$SLURM_ARRAY_TASK_ID \
		RGLB=TruSeq_stranded_mRNA RGPL=illumina \
		RGPU=NextSeq RGSM=$SLURM_ARRAY_TASK_ID

# Picard Markduplicates (for amplicon analyses, this step may be skipped)
java -jar $EBROOTPICARD/picard.jar MarkDuplicates I=$OutFOLDER/$InBAM.rg.sort.bam O=$OutFOLDER/$InBAM.rg.sort.dedupped.bam CREATE_INDEX=true \
			VALIDATION_STRINGENCY=SILENT M=$OutFOLDER/$InBAM.duplicate.metrics.txt 

# index bam files
samtools index $OutFOLDER/$InBAM.rg.sort.dedupped.bam

