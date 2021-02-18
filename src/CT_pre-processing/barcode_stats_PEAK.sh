#!/bin/bash
#SBATCH -A snic2020-5-572
#SBATCH -n 4
#SBATCH -p core
#SBATCH -t 1-0

#####
#####



ml bioinfo-tools samtools BEDTools

BASE_DIR="/proj/uppstore2017150/private/mukund_playground/ATAC/10x_data"

bedtools intersect -abam ${BASE_DIR}/$1/outs/possorted_bam.bam -b $1/$1_SEACR.relaxed.bed | \
	samtools view -f2 | \
		awk -f scripts/get_cell_barcode.awk | \
			sed 's/CB:Z://g' | sort | uniq -c > barcode_statistics/$1_relaxed.peak_barcodes.txt && [[ -s barcode_statistics/$1_relaxed.peak_barcodes.txt ]]




bedtools intersect -abam ${BASE_DIR}/$1/outs/possorted_bam.bam -b $1/$1_SEACR.stringent.bed | \
	samtools view -f2 | \
		awk -f scripts/get_cell_barcode.awk | \
			sed 's/CB:Z://g' | sort | uniq -c > barcode_statistics/$1_stringent.peak_barcodes.txt && [[ -s barcode_statistics/$1_stringent.peak_barcodes.txt ]]






