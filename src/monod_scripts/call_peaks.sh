#!/bin/bash

/home/mukund/miniconda3/envs/episcanpy/bin/python src/split_fragments_by_cluster.py \
	--fragments /proj/tmp/mukund/ct/P22551_1019_v2/outs/frags.tsv.gz \
	--barcodes ../211103_scCT_barcodes_by_leiden2.0_cluster_LSI-1.csv \
	--output data/cluster_fragments/LSI_1_scCT_fragments_ 

cd data/cluster_fragments

# call peaks with ATAC-compatible parameters
for FILE in *.bed
do
	~/miniconda3/envs/episcanpy/bin/macs2 callpeak \
	-t $FILE \
	-f BED \
	-g hs \
	-n $FILE \
	-q 0.05 \
	--shift -100 \
	--extsize 200 \
	--nomodel \
	--call-summits --keep-dup=1 \
	--outdir ../peaks
	
done


#########################
#### ONLY FOR LSI-1 #####
#########################

cd ../peaks

### resize peaks to 150bp centred on the summit
for FILE in *summits.bed
do
	/data/bin/bedtools2/bin/bedtools slop -i $FILE -g ../../ref/hg38.chrom.sizes -b 75 | \
	awk 'OFS = "\t" {print $1,$2,$3,NR}' > PADDED_${FILE}
done


## merge to get final set of peaks (use to build the feature matrix)
cat PADDED*bed | sort -k1,1 -k2,2n | \
/data/bin/bedtools2/bin/bedtools merge -i stdin > 211103_scCT_merged_leiden2.0_cluster_peaks.bed

cd ../..

PYTHONPATH="/home/mukund/miniconda3/envs/episcanpy/lib/python3.6"
PYTHONHOME="/home/mukund/miniconda3/envs/episcanpy"

/home/mukund/miniconda3/envs/episcanpy/bin/python src/build_custom_feature_matrix.py \
	-f /proj/tmp/mukund/ct/P22551_1019_v2/outs/frags.tsv.gz \
	-p data/peaks/211103_scCT_merged_leiden2.0_cluster_peaks.bed \
	-c /proj/tmp/mukund/ct/P22551_1019_v2/outs/singlecell_EPI.csv \
	-o ../211103_scCT_mergedPeak_matrix.h5ad #\
#	-m ref/sample_metadata_COMPLETE.csv \
#	-l data/feature_matrices/211101_3x1K_merged_peak_matrix_LOADED.h5ad
