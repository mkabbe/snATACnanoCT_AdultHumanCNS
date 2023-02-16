#!/bin/bash

### Export BW tracks ###


cd /data/proj/GCB_MK/CZI/episcanpy_analysis/AGG_ATAC_210218

## Split fragments by celltype
/home/mukund/miniconda3/envs/episcanpy/bin/python src/split_fragments_by_cluster.py \
       --fragments /data/proj/GCB_MK/10XATAC/data/scatac/MAIN_AGGREGATE_ALL_CELLS/fragments.tsv.gz \
       --barcodes data/221210_annot_barcodes_hiQ.csv \
       --output data/cluster_fragments/CELLTYPE_221210/

#cd data/cluster_fragments/OL_SUB/
#cd data/cluster_fragments/OL_TISSUE/
#cd data/cluster_fragments/OL_SEX/
#cd data/cluster_fragments/OL_AGE/
#cd data/cluster_fragments/ALL_CELL_TYPES/
#cd data/cluster_fragments/OL_CCA_TISSUE/
#cd data/cluster_fragments/OL_OLG_TISSUE/

#cd data/cluster_fragments/OL_typeTissue/
#cd data/cluster_fragments/OL_CASENO/
#cd data/cluster_fragments/OPC_RegionAgeCat/
#cd data/cluster_fragments/OL_RegionAgeCat/

#cd data/cluster_fragments/nonOL_region-cells/

cd data/cluster_fragments/CELLTYPE_221210/


## Convert to BAM and make BigWig file

mkdir -p BAM_FILES
mkdir -p BW_FILES

for FILE in *bed
do
         bedtools bedtobam -i $FILE -g ../../../ref/hg38.chrom.sizes > BAM_FILES/$FILE.bam
         ~/miniconda3/envs/episcanpy/bin/samtools index BAM_FILES/$FILE.bam
#         ~/miniconda3/envs/episcanpy/bin/bamCoverage -b BAM_FILES/$FILE.bam -o BW_FILES/$FILE.rpkm.bw -p max --normalizeUsing RPKM
done

#computeMatrix reference-point -R hox_regions/hoxa_promoters.bed -S BW_FILES/*bw -a 2000 -b 2000 -o hox_regions/hox_coverage_Age.matrix
#plotProfile -m hox_regions/hox_coverage_Age.matrix -o hox_regions/plots/hox_coverage_age.pdf -T HOXA_coverage --perGroup --plotHeight 20 --plotWidth 30
