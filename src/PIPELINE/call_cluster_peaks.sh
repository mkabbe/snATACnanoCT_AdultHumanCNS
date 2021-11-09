
## ASSUMES YOU ARE IN THE PIPELINE DIRECTORY

# has 1 command line argument

cd data/cluster_fragments

## CLA
MERGED_PEAK_BED=$1


## call peaks with ATAC compatible parameters
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

## resize peaks to 150bp centred on the summit
for FILE in *summits.bed
do
	/data/bin/bedtools2/bin/bedtools slop -i $FILE -g ../../ref/hg38.chrom.sizes -b 75 | \
	awk 'OFS = "\t" {print $1,$2,$3,NR}' > PADDED_${FILE}
done


## merge to get final set of peaks (use to build the feature matrix)
cat PADDED*bed | sort -k1,1 -k2,2n | \
/data/bin/bedtools2/bin/bedtools merge -i stdin > MERGED_PEAK_BED
