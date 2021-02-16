#!/bin/bash
#SBATCH -A snic2017-7-345
#SBATCH -n 4
#SBATCH -p core
#SBATCH -t 1-0

ml bioinfo-tools samtools BEDTools

#samtools view possorted_bam.bam | grep CB:Z: | sed 's/.*CB:Z:\([ACGT]*\).*/\1/' | sort | uniq -c > reads_per_barcode_on_core


bedtools intersect -abam CT_20_MK01_k4me3/possorted_bam.bam -b CT_20_MK01_k4me3/CT_20_MK01_k4me3_peaks.narrowPeak -u | \
samtools view -f2 | \
awk -f scripts/get_cell_barcode.awk | \
sed 's/CB:Z://g' | sort | uniq -c > CT_20_MK01_k4me3/barcode_statistics/peaks_barcodes_narrow.txt && [[ -s CT_20_MK01_k4me3/barcode_statistics/peaks_barcodes_narrow.txt ]]







#rule barcode_statistics_peaks:
#    input:
#        bam          = lambda wildcards: config['samples'][wildcards.sample]['cellranger_out'] + '/outs/possorted_bam.bam',
#        peaks_broad  = "results/{sample}/macs/broad/{sample}_peaks.broadPeak",
#        peaks_narrow = "results/{sample}/macs/narrow/{sample}_peaks.narrowPeak"
#    output:
#        narrow = "results/{sample}/barcode_statistics/peaks_barcodes_narrow.txt",
#        broad  = "results/{sample}/barcode_statistics/peaks_barcodes_broad.txt"
#    params:
#        scripts    = os.path.dirname(workflow.basedir) + '/scripts',
#    shell:
#      #  "set +o pipefail; "
#        "bedtools intersect -abam {input.bam} -b {input.peaks_broad} -u | samtools view -f2 | "
#        " awk -f {params.scripts}/get_cell_barcode.awk | sed 's/CB:Z://g' | sort | uniq -c > {output.broad} && [[ -s {output.broad} ]] ; "
#        " bedtools intersect -abam {input.bam} -b {input.peaks_narrow} -u | samtools view -f2 | "
#        " awk -f {params.scripts}/get_cell_barcode.awk | sed 's/CB:Z://g' | sort | uniq -c > {output.narrow} && [[ -s {output.narrow} ]] ;"


