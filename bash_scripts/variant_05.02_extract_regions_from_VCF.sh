#!/usr/bin/bash

OUTDIR=variant_05.02_extract_regions_from_VCF

mkdir -p $OUTDIR

SEEDFILE=../various_files/regions_for_SNPs_of_interest.txt
SEED_STR=$(sed -n ${SLURM_ARRAY_TASK_ID}p $SEEDFILE)
IFS=' ' read -r -a SEED_ARRAY <<< "$SEED_STR"

GENE="${SEED_ARRAY[0]}"
GENENAME="${SEED_ARRAY[1]}"
CHR="${SEED_ARRAY[2]}"
START="${SEED_ARRAY[3]}"
STOP="${SEED_ARRAY[4]}"

INFILE=variant_05_GenotypeGVCFs.snps.filtered.vcf.gz

OUTFILE=$OUTDIR/${GENE}_${GENENAME}.${CHR}:${START}-${STOP}.vcf
.vcf

ml VCFtools
time vcftools --gzvcf $INFILE --chr $CHR --from-bp $START --to-bp $STOP --recode --recode-INFO-all --out $OUTFILE
