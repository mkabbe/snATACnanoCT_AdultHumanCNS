#!/bin/bash
#SBATCH -A snic2017-7-345
#SBATCH -n 4
#SBATCH -p core
#SBATCH -t 1-0

ml bioinfo-tools samtools BEDTools

#samtools view possorted_bam.bam | grep CB:Z: | sed 's/.*CB:Z:\([ACGT]*\).*/\1/' | sort | uniq -c > reads_per_barcode_on_core

samtools view -f2 CT_20_MK01_k4me3/possorted_bam.bam | \
awk -f scripts/get_cell_barcode.awk | \
sed 's/CB:Z://g' | sort | uniq -c > CT_20_MK01_k4me3/barcode_statistics/all_barcodes.txt  && [[ -s CT_20_MK01_k4me3/barcode_statistics/all_barcodes.txt ]]





#rule barcode_statistics_all:
#  input:
#     bam       = lambda wildcards: config['samples'][wildcards.sample]['cellranger_out'] + '/outs/possorted_bam.bam',
#  output:
#    all_bcd    = "results/{sample}/barcode_statistics/all_barcodes.txt"
#  params:
#    scripts    = os.path.dirname(workflow.basedir) + '/scripts',
#  shell:
#    #" set +o pipefail; "
#    " samtools view -f2 {input.bam}| "
#    " awk -f {params.scripts}/get_cell_barcode.awk | sed 's/CB:Z://g' | sort | uniq -c > {output.all_bcd} && [[ -s {output.all_bcd} ]] "


