#!/bin/bash
#SBATCH -A snic2020-5-572
#SBATCH -n 4
#SBATCH -p core
#SBATCH -t 1-0

ml bioinfo-tools samtools BEDTools

#samtools view possorted_bam.bam | grep CB:Z: | sed 's/.*CB:Z:\([ACGT]*\).*/\1/' | sort | uniq -c > reads_per_barcode_on_core

samtools view -f2 CT_20_MK01_k4me3/possorted_bam.bam | \
awk -f scripts/get_cell_barcode.awk | \
sed 's/CB:Z://g' | sort | uniq -c > CT_20_MK01_k4me3/barcode_statistics/all_barcodes.txt  && [[ -s CT_20_MK01_k4me3/barcode_statistics/all_barcodes.txt ]]


################
################

ml bioinfo-tools samtools BEDTools

BASE_DIR="/proj/uppstore2017150/private/mukund_playground/ATAC/10x_data"

samtools view -f2 ${BASE_DIR}/$1/outs/possorted_bam.bam | \
	awk -f scripts/get_cell_barcode.awk | \
		sed 's/CB:Z://g' | sort | uniq -c > barcode_statistics/$1_all_barcodes.txt && [[ -s barcode_statistics/$1_all_barcodes.txt ]]

