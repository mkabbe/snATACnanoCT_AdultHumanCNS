#!/bin/bash


## Move into the project directory before running this script!
#       e.g. Be in P28305/


# Change project number and # samples in for loop as needed


for sample in 1001 1002 1003 1004 1005 1006
do

cd P28305_${sample}

mv 02*FASTQ/*/*/* .

echo $sample
prefix=$(echo $sample | cut -d0 -f3)
cat *_I1_* > P28305_${sample}_S${prefix}_I1_001.fastq.gz
cat *_R1_* > P28305_${sample}_S${prefix}_R1_001.fastq.gz
cat *_R2_* > P28305_${sample}_S${prefix}_R2_001.fastq.gz
cat *_R3_* > P28305_${sample}_S${prefix}_R3_001.fastq.gz

rm *L00*

cd ..