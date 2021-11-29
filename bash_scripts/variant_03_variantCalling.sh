#!/usr/bin/bash

ml GATK

SAMPLE=`sed -n ${SLURM_ARRAY_TASK_ID}p  sampleNames.txt`

GENOME_FASTA=Anopheles-coluzzii-Ngousso_CHROMOSOMES_AcolN1.fa

InFOLDER="variant_02_bamPreProcessing"
OutFOLDER="variant_03_variantCalling"

InBAM=$SAMPLE.merged.Aligned.sortedByCoord.out.bam.rg.sort.dedupped.bam

# GATK SplitNCigarReads
# At this step we also add one important tweak: we need to reassign mapping qualities, 
# because STAR assigns good alignments a MAPQ of 255 (which technically means "unknown" 
# and is therefore meaningless to GATK). So we use the GATK's ReassignOneMappingQuality read filter 
# to reassign all good alignments to the default value of 60. This is not ideal, and we hope that in 
# the future RNAseq mappers will emit meaningful quality scores, but in the meantime this is the best we can do. 
# In practice we do this by adding the ReassignOneMappingQuality read filter to the splitter command.

# gatk --java-options "-Xmx2G" SplitNCigarReads -R $GENOME_FASTA -I $InFOLDER/$InBAM -O $OutFOLDER/$InBAM.split.bam 

# GATK HaplotypeCaller 
gatk HaplotypeCaller --java-options "-Xmx2G" -R $GENOME_FASTA -I $OutFOLDER/$InBAM.split.bam \
    -stand-call-conf 20 -O $OutFOLDER/$InBAM.split.bam.g.vcf -ERC GVCF

