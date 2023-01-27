#!/bin/bash
#SBATCH -A snic2021-22-1020
#SBATCH -n 4
#SBATCH -p core
#SBATCH -t 1-0
##SBATCH -p devcore
##SBATCH -t 60

ml bioinfo-tools samtools BEDTools

NAME=$1  ## e.g. P27308_1001
ANTIBODY=$2  ## e.g. H3K27ac
WIDTH=$3  ## narrow (or broad)

mkdir -p BAM/BCD_STATS_{ALL,PEAK}

## All fragments
samtools view -f2 BAM/"$NAME"_"$ANTIBODY"/possorted_bam.bam | awk -f src/get_cell_barcode.awk | \
sed 's/CB:Z://g' | sort | uniq -c > BAM/BCD_STATS_ALL/"$NAME"_"$ANTIBODY"_all_barcodes.txt && [[ -s BAM/BCD_STATS_ALL/"$NAME"_"$ANTIBODY"_all_barcodes.txt  ]]


## Fragments in peaks
bedtools intersect -abam BAM/"$NAME"_"$ANTIBODY"/possorted_bam.bam -b BAM/PEAKS/"$NAME"_"$ANTIBODY"_peaks."$WIDTH"Peak -u | \
samtools view -f2 | awk -f src/get_cell_barcode.awk | \
sed 's/CB:Z://g' | sort | uniq -c > BAM/BCD_STATS_PEAK/"$NAME"_"$ANTIBODY"_"$WIDTH"peak_barcodes.txt && [[ -s BAM/BCD_STATS_PEAK/"$NAME"_"$ANTIBODY"_"$WIDTH"peak_barcodes.txt ]]


