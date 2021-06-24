#!/bin/bash

### Export BW tracks ###


cd /data/proj/GCB_MK/CZI/episcanpy_analysis/AGG_ATAC_210218

## Split fragments by celltype
/home/mukund/miniconda3/envs/episcanpy/bin/python src/split_fragments_by_cluster.py \
       --fragments /proj/tmp/mukund/atac/AGG_ATAC_210507/outs/fragments.tsv.gz \
       --barcodes data/210622_cellType_barcodes_LSI-2.csv \
       --output data/cluster_fragments/all_

cd data/cluster_fragments


## Convert to BAM and make BigWig file

for FILE in *bed
do
        bedtools bedtobam -i $FILE -g ../../ref/hg38.chrom.sizes > BAM_FILES/$FILE.bam
        samtools index BAM_FILES/$FILE.bam
        bamCoverage -b BAM_FILES/$FILE.bam -o BW_FILES/$FILE.bw -p max
done