#!/usr/bin/bash

ml GATK

gatk --java-options "-Xmx4g" GenotypeGVCFs \
   -R Anopheles-coluzzii-Ngousso_CHROMOSOMES_AcolN1.fa \
   -V gendb://variant_04_GenomicsDBImport \
   -O variant_05_GenotypeGVCFs.vcf.gz

