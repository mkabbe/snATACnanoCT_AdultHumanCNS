
### Do Iterative clustering on cluster-specific peaks


## specify cluster name for running on condor separately
CLUSTER=$1
mkdir -p __tmp__$CLUSTER ## temp folders for each cluster
cd __tmp__$CLUSTER

## iterate through all chromosomes
CHROMS=('1' '2' '3' '4' '5' '6' '7' '8' '9' '10' '11' '12' '13' '14' '15' '16' '17' '18' '19' '20' '21' '22' 'X' 'Y')

for i in "${CHROMS[@]}"
do
	echo "** chr$i **" >> iter.log 
	SECONDS=0 
	
	## Select chromosome, sort by summit start position and then by score
	grep chr$i ../data/peaks/lsi1_leiden_$CLUSTER.bed_summits.bed | sort -k1,1V -k2,2n -k3,5n > TMP.bed
	
	total_peaks=$(cat TMP.bed | wc -l)	
	echo "** chr$i : $total_peaks total peaks **" >> iter.log	
	
	## Expand summit to 500bp
	bedtools slop -i TMP.bed -g ../ref/hg38.chrom.sizes -b 250 > PEAK.bed

	## Iterative removal of less significant peaks
	while [[ $(cat PEAK.bed | wc -l) -ne 0 ]];
	do
		head -n1  PEAK.bed > top.bed ## get TOP peak in file
		cat top.bed >> IterPeak.bed ## store in final peak file
		bedtools intersect -v -sorted  -a PEAK.bed -b top.bed > TMP.bed ## find peaks that DO NOT overlap with top peak
		mv TMP.bed PEAK.bed ## replace peak file

		## Track number of peaks left
		peaks_left=$(cat PEAK.bed | wc -l)
		if [[ $(($peaks_left % 5000)) -eq 0 ]]; ## print every 5000 peaks
		then
			echo "$peaks_left peaks left" >> iter.log
		fi
	done
	
	## print total numer of filtered peaks remaining
	filtered_peaks=$(cat IterPeak.bed | wc -l) 
	echo "** chr$i : $filtered_peaks / $total_peaks peaks **" >> iter.log
	ELAPSED="$(($SECONDS / 3600))hrs $((($SECONDS / 60) % 60))min $(($SECONDS % 60))sec"
	echo "** chr$i :  $ELAPSED **" >> iter.log
	
	## add to the final peak file
	cat IterPeak.bed >> ../data/peaks/FILTERED_lsi1_leiden_$CLUSTER.bed_summits.bed
done

#cd ..
#rm -r __tmp__$CLUSTER #delete tmp folder
