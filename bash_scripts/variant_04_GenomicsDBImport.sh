#!/usr/bin/bash

ml GATK

# request 6Gb
# MPORTANT: The -Xmx value the tool is run with should be less than the total amount of physical memory available by at least a few GB, as the native TileDB library requires additional memory on top of the Java memory. Failure to leave enough memory for the native code can result in confusing error messages!

# https://gatk.broadinstitute.org/hc/en-us/articles/360035889971
TILEDB_DISABLE_FILE_LOCKING=1

# The --genomicsdb-workspace-path must point to a non-existent or empty directory.
outDIR=variant_04_GenomicsDBImport
rm -rf $outDIR
#mkdir -p $outDIR 

ls variant_03_variantCalling/*.g.vcf |  awk 'BEGIN{FS="."}{ print $1"\t"$0 }' | sed 's|^variant_03_variantCalling/||g' > sample_map.tmp.txt
grep "##contig" `head -n 1 sample_map.tmp.txt | cut -f 2` |  cut -d , -f 1 | cut -d = -f 3 > chr_list.tmp.list

gatk --java-options "-Xmx4g -Xms4g" \
       GenomicsDBImport \
       --genomicsdb-workspace-path $outDIR \
       --batch-size 50 \
       -L chr_list.tmp.list \
       --sample-name-map sample_map.tmp.txt \
       --tmp-dir=$TMPDIR \
       --reader-threads 1

rm chr_list.tmp.list
rm sample_map.tmp.txt
