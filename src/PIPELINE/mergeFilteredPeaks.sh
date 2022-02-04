#!/bin/bash -l
#SBATCH -A snic2021-22-1020
#SBATCH -p core
#SBATCH -n 16
#SBATCH -t 2-0
#SBATCH -J IterativePeaks_CHROM_NAME
#SBATCH -e logs/err_CHROM_NAME.err -o logs/out_CHROM_NAME.out
#SBATCH --mail-user=mukund.kabbe@ki.se --mail-type=BEGIN,END

module load bioinfo-tools BEDTools/2.29.2 ## rackham specific

### Do Iterative clustering on cluster-specific peaks

## specify cluster name for running on condor separately
CHROM=$1
mkdir -p __tmp__$CHROM ## temp folders for each cluster
mkdir -p logs/Iter
cd __tmp__$CHROM

SECONDS=0 

## Select chromosome, sort by summit start position and then by score
#grep -w chr$i ../data/peaks/lsi1_leiden_$CLUSTER.bed_summits.bed | sort -k1,1V -k2,2n -k3,5n > TMP.bed

cp ../data/peaks/all_FILTERED_chr${CHROM}_peaks.bed PEAK.bed

total_peaks=$(cat PEAK.bed | wc -l)	
echo "** chr$CHROM : $total_peaks total peaks **" >> iter.log	

## Expand summit to 500bp
#bedtools slop -i TMP.bed -g ../ref/hg38.chrom.sizes -b 250 > PEAK.bed

## Iterative removal of less significant peaks
while [[ $(cat PEAK.bed | wc -l) -ne 0 ]];
do
	head -n1  PEAK.bed > top.bed ## get TOP peak in file
	cat top.bed >> IterPeak.bed ## store in final peak file
	bedtools intersect -v -sorted  -a PEAK.bed -b top.bed > TMP.bed ## find peaks that DO NOT overlap with top peak
	mv TMP.bed PEAK.bed ## replace peak file

	## Track number of peaks left
	peaks_left=$(cat PEAK.bed | wc -l)
	if [[ $(($peaks_left % 50000)) -eq 0 ]]; ## print every 50000 peaks
	then
		echo "$peaks_left peaks left" >> iter.log
	fi
done
	
## print total numer of filtered peaks remaining
filtered_peaks=$(cat IterPeak.bed | wc -l) 
echo "** chr$CHROM : $filtered_peaks / $total_peaks peaks **" >> iter.log
ELAPSED="$(($SECONDS / 3600))hrs $((($SECONDS / 60) % 60))min $(($SECONDS % 60))sec"
echo "** chr$CHROM :  $ELAPSED **" >> iter.log

mv IterPeak.bed ../data/peaks/REFILTERED_all_chr${CHROM}_peaks.bed

## add to the final peak file
#cat IterPeak.bed >> ../data/peaks/FILTERED_lsi1_leiden_$CLUSTER.bed_summits.bed
#rm IterPeak.bed ## clear peak file for next chromosome

#mv iter.log ../logs/Iter/iter_log_$CHROM.log ## store logs in different folder
#cd ..
#rm -r __tmp__$CHROM #delete tmp folder
