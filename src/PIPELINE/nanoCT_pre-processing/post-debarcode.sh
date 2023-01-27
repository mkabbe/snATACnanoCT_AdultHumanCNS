#!/bin/bash


## This moves and renames the files that debarcode.py outputs.
## Run before running cellranger-atac. 

## Change project number and # samples in for loop as needed

## TATAGCCT - H3K27ac
## ATAGAGGC - H3K27me3


cd /data/proj/GCB_MK/scCT/nanoCT_EAE/fastq/debarcoded


for sample in 1001 1002 1003 1004 1005 1006
do

        mkdir split/P28305_${sample}_H3K27me3
        mkdir split/P28305_${sample}_H3K27ac

        mv P28305_${sample}/barcode_TATAGCCT/*R1* split/P28305_${sample}_H3K27ac/P28305_${sample}_H3K27ac_S1_R1_001.fastq.gz
        mv P28305_${sample}/barcode_TATAGCCT/*R2* split/P28305_${sample}_H3K27ac/P28305_${sample}_H3K27ac_S1_R2_001.fastq.gz
        mv P28305_${sample}/barcode_TATAGCCT/*R3* split/P28305_${sample}_H3K27ac/P28305_${sample}_H3K27ac_S1_R3_001.fastq.gz

        mv P28305_${sample}/barcode_ATAGAGGC/*R1* split/P28305_${sample}_H3K27me3/P28305_${sample}_H3K27me3_S1_R1_001.fastq.gz
        mv P28305_${sample}/barcode_ATAGAGGC/*R2* split/P28305_${sample}_H3K27me3/P28305_${sample}_H3K27me3_S1_R2_001.fastq.gz
        mv P28305_${sample}/barcode_ATAGAGGC/*R3* split/P28305_${sample}_H3K27me3/P28305_${sample}_H3K27me3_S1_R3_001.fastq.gz

done