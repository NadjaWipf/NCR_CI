#!/usr/bin/bash

INDIR=variant_05.03_classify_regions_SNPs
OUTDIR=variant_05.04_regions_VCF_table

mkdir -p $OUTDIR

SEEDFILE=../various_files/regions_for_SNPs_of_interest.txt
SEED_STR=$(sed -n ${SLURM_ARRAY_TASK_ID}p $SEEDFILE)
IFS=' ' read -r -a SEED_ARRAY <<< "$SEED_STR"

GENE="${SEED_ARRAY[0]}"
GENENAME="${SEED_ARRAY[1]}"
CHR="${SEED_ARRAY[2]}"
START="${SEED_ARRAY[3]}"
STOP="${SEED_ARRAY[4]}"

INPUTVCF=$INDIR/${GENE}_${GENENAME}.${CHR}:${START}-${STOP}.vcf.recode.vcf

OUTPUTVCF=$OUTDIR/${GENE}_${GENENAME}.${CHR}:${START}-${STOP}.filtered.vcf
OUTPUTTABLE=$OUTDIR/${GENE}_${GENENAME}.${CHR}:${START}-${STOP}.filtered.tsv

ml GATK

gatk VariantFiltration \
-V $INPUTVCF \
-O $OUTPUTVCF \
--set-filtered-genotype-to-no-call \
--genotype-filter-expression "DP < 3" \
--genotype-filter-name "DP-3" \
--genotype-filter-expression "GQ < 10" \
--genotype-filter-name "GQ-10"

gatk VariantsToTable \
     -V $OUTPUTVCF \
     -F CHROM -F POS -F ID -F REF -F ALT -F QUAL -GF GT \
     -O $OUTPUTTABLE

