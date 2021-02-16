#!/bin/bash
#SBATCH -A snic2020-5-572
#SBATCH -p devcore
#SBATCH -n 4
#SBATCH -t 60 
#### include next time --mem=64G

module load bioinfo-tools samtools BEDTools

samtools sort -n $1_possorted_bam.bam > $1_namesorted_bam.bam
bedtools bamtobed -bedpe -i $1_namesorted_bam.bam > $1.bed
awk '$1==$4 && $6-$2 < 1000 {print $0}' $1.bed > $1.clean.bed
cut -f 1,2,6 $1.clean.bed | sort -k1,1 -k2,2n -k3,3n > $1.fragments.bed
awk 'FNR==NR {a[$1]; next} $1 in a' hg38.chrom.sizes $1.fragments.bed | bedtools genomecov -bg -i stdin -g hg38.chrom.sizes > $1.fragments.bedgraph

bash SEACR_1.3.sh $1.fragments.bedgraph 0.01 non stringent $1_SEACR

