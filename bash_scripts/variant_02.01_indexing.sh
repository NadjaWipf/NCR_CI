#!/usr/bin/bash

GENOME_FASTA=Anopheles-coluzzii-Ngousso_CHROMOSOMES_AcolN1.fa

ml picard
ml SAMtools
java -jar $EBROOTPICARD/picard.jar CreateSequenceDictionary R=$GENOME_FASTA O=$GENOME_FASTA.dict
samtools faidx $GENOME_FASTA

