#!/bin/bash
#SBATCH -A snic2020-5-572
#SBATCH -n 4
#SBATCH -p core
#SBATCH -t 1-0

################
################

ml bioinfo-tools samtools BEDTools

BASE_DIR="/proj/uppstore2017150/private/mukund_playground/ATAC/10x_data"

samtools view -f2 ${BASE_DIR}/$1/outs/possorted_bam.bam | \
	awk -f scripts/get_cell_barcode.awk | \
		sed 's/CB:Z://g' | sort | uniq -c > barcode_statistics/$1_all_barcodes.txt && [[ -s barcode_statistics/$1_all_barcodes.txt ]]

