#!/bin/bash

cd /data/proj/GCB_MK/CZI/episcanpy_analysis/AGG_ATAC_210218

## Split fragments by celltype
#/home/mukund/miniconda3/envs/episcanpy/bin/python src/split_fragments_by_cluster.py \
#       --fragments /data/proj/GCB_MK/10XATAC/data/scatac/MAIN_AGGREGATE_ALL_CELLS/fragments.tsv.gz \
#       --barcodes data/221210_annot_barcodes_hiQ.csv \
#       --output data/cluster_fragments/CELLTYPE_221210/


cd data/cluster_fragments/CELLTYPE_221210/
mkdir -p BAM_FILES
mkdir -p BW_FILES
mkdir -p PEAKS


# call peaks with ATAC-compatible parameters
for FILE in *.bed
do
	~/miniconda3/envs/episcanpy/bin/macs2 callpeak \
	-t $FILE \
	-f BED \
	-g hs \
	-n $FILE \
	-q 0.01 \
	--shift -100 \
	--extsize 200 \
	--nomodel \
	--call-summits --keep-dup=1 \
	--outdir PEAKS	
done


## Convert to BAM
for FILE in *bed
do
         bedtools bedtobam -i $FILE -g ../../../ref/hg38.chrom.sizes > BAM_FILES/$FILE.bam
         ~/miniconda3/envs/episcanpy/bin/samtools index BAM_FILES/$FILE.bam
done


