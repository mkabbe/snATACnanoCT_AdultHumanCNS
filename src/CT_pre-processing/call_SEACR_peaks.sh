#!/bin/bash
#SBATCH -A snic2020-5-572
#SBATCH -p core
#SBATCH -n 4
#SBATCH -t 1-0 
#SBATCH --mem=64G


module load bioinfo-tools samtools BEDTools

BASE_DIR="/proj/uppstore2017150/private/mukund_playground/ATAC/10x_data"


samtools sort -n ${BASE_DIR}/$1/outs/possorted_bam.bam > $1/$1_namesorted_bam.bam 
bedtools bamtobed -bedpe -i $1/$1_namesorted_bam.bam > $1/$1.bed
awk '$1==$4 && $6-$2 < 1000 {print $0}' $1/$1.bed > $1/$1.clean.bed
cut -f 1,2,6 $1/$1.clean.bed | sort -k1,1 -k2,2n -k3,3n > $1/$1.fragments.bed
awk 'FNR==NR {a[$1]; next} $1 in a' hg38.chrom.sizes $1/$1.fragments.bed | bedtools genomecov -bg -i stdin -g hg38.chrom.sizes > $1/$1.fragments.bedgraph

bash scripts/SEACR_1.3.sh $1/$1.fragments.bedgraph 0.1 non relaxed $1_SEACR && mv $1_SEACR.relaxed.bed $1/
bash scripts/SEACR_1.3.sh $1/$1.fragments.bedgraph 0.1 non stringent $1_SEACR && mv $1_SEACR.stringent.bed $1/
