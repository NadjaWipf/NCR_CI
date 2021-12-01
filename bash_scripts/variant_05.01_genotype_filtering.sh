#!/usr/bin/bash

ml GATK

# let's restrict ourselves to snps
gatk --java-options "-Xmx2g" SelectVariants \
    -V variant_05_GenotypeGVCFs.vcf.gz \
    -select-type SNP \
    -O variant_05_GenotypeGVCFs.snps.vcf.gz

gatk --java-options "-Xmx2g" VariantFiltration \
    -V variant_05_GenotypeGVCFs.snps.vcf.gz \
    -filter "QUAL < 30.0" --filter-name "QUAL30" \
    -filter "MQ < 40.0" --filter-name "MQ40" \
    -filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
    -filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
    -O variant_05_GenotypeGVCFs.snps.filtered.vcf.gz

