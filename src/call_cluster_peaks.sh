#!/bin/bash 

## if running on cluster:
# ml bioinfo-tools samtools MACS/2.1.0 BEDTools

# call peaks with ATAC-compatible parameters
for FILE in data/epi/kmeans_fragments/*.bed
do
        macs2 callpeak \
        -t $FILE \
        -f BED \
        -g hs \
        -n $FILE \
        -q 0.05 \
        --shift -100 \
        --extsize 200 \
        --nomodel \
        --call-summits --keep-dup=1 \
        --outdir data/epi/cluster_peaks
		
done

## resize peaks to 150bp centred on the summit
for FILE in *summits.bed
do
       bedtools slop -i $FILE -g ../../hg38.chrom.sizes -b 75 | \
	   awk 'OFS = "\t" {print $1,$2,$3,NR}' > data/epi/cluster_peaks/PADDED_${FILE}
done

## merge to get final set of peaks (use to build the feature matrix)
cat data/epi/kmeans_fragments/PADDED*bed | sort -k1,1 -k2,2n | \
bedtools merge -i stdin > data/epi/kmeans_fragments/merged_kmeans_cluster_peaks.bed