#!/bin/bash
#SBATCH -A snic2021-22-1020
#SBATCH -p core
#SBATCH -n 12
#SBATCH -t 1-0
#SBATCH --mem=16G

ml bioinfo-tools MACS/2.2.6

NAME=$1     ## e.g. P27308_1001
ANTIBODY=$2     ## e.g. H3K27ac 

mkdir -p BAM/PEAKS

macs2 callpeak -t BAM/"$NAME"_"$ANTIBODY"/possorted_bam.bam -g hs -f BAMPE -n "$NAME"_"$ANTIBODY" --outdir BAM/PEAKS --keep-dup=1 --llocal 100000 --min-length 1000 --max-gap 1000 2>&1



##broad - H3K27me3
#macs2 callpeak -t BAM/P27308_1001_H3K27me3/possorted_bam.bam -g mm -f BAMPE -n P27308_1001_H3K27me3 --outdir BAM/PEAKS/macs_broad/ --keep-dup=1 --llocal 100000 --min-length 1000 --max-gap 1000 --broad-cutoff=0.1 --broad 2>&1

