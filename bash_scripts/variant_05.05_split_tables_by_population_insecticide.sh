#!/usr/bin/bash

ls variant_05.04_regions_VCF_table/*.filtered.tsv > tmp
for IN in `cat tmp`
do
cut -f 1,2,3,4,5,6,7,58,59,60,61 $IN >$IN.Agb_M.tsv
cut -f 1,2,3,4,5,6,17,28,39,50,57 $IN >$IN.Agb_D.tsv
cut -f 1,2,3,4,5,6,8,9,10,11,12 $IN >$IN.Agb_C.tsv
cut -f 1,2,3,4,5,6,41,42,43,44,45 $IN >$IN.Tia_C.tsv
cut -f 1,2,3,4,5,6,35,36,37,38,40 $IN >$IN.Tia_M.tsv
cut -f 1,2,3,4,5,6,30,31,32,33,34 $IN >$IN.Tia_D.tsv
cut -f 1,2,3,4,5,6,24,25,26,27,29 $IN >$IN.Dab_C.tsv
cut -f 1,2,3,4,5,6,19,20,21,22,23 $IN >$IN.Dab_M.tsv
cut -f 1,2,3,4,5,6,13,14,15,16,18 $IN >$IN.Dab_D.tsv
cut -f 1,2,3,4,5,6,46,47,48,49,51 $IN >$IN.Ng_C.tsv
cut -f 1,2,3,4,5,6,52,53,54,55,56 $IN >$IN.Ma_C.tsv
done

