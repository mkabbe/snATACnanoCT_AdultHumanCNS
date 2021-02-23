#!/bin/bash 


# run this on the output of the python notebook before running ArchR


#sed '/^[[:space:]]*$/d' ../CT_fragments/MK02_H3K27ac_FILTERED.fragments.tsv > tmp && mv tmp ../CT_fragments/MK02_H3K27ac_FILTERED.fragments.tsv
#bgzip < ../CT_fragments/MK02_H3K27ac_FILTERED.fragments.tsv > ../CT_fragments/MK02_H3K27ac_FILTERED.fragments.tsv.gz

sed '/^[[:space:]]*$/d' data/CT_processed_fragments/P20057_1004_FILTERED.fragments.tsv > tmp 
mv tmp data/CT_processed_fragments/P20057_1004_FILTERED.fragments.tsv
bgzip < data/CT_processed_fragments/P20057_1004_FILTERED.fragments.tsv > data/CT_processed_fragments/P20057_1004_FILTERED.fragments.tsv.gz

