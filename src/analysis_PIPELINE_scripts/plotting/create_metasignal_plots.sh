#!/bin/bash


conda activate pycistopic

computeMatrix reference-point \
	-R sorted_OLG_H3K27ac_fragment.bed_peaks.broadPeak \
		-S BW_FILES/sorted_OLG_H3K27*.bw \
			-o OLG_H3K27ac_peaks.matrix \
				-b 10000 -a 10000



computeMatrix reference-point -R sorted_OLG_H3K27ac_fragment.bed_peaks.narrowPeak -S BW_FILES/sorted_OLG_H3K27*.bw -o OLG_H3K27ac_peaks.matrix -b 10000 -a 10000

## @ H3K27me3 peaks

# H3K27me3 signal
computeMatrix reference-point -R H3K27me3.peaks -S BW_FILES/sorted_OLG_H3K27ac*bw -o OLG_H3K27ac_signal_me3_Peaks.matrix.gz -b 10000 -a 10000 -p 32
plotHeatmap -m OLG_H3K27ac_signal_me3_Peaks.matrix.gz -o OLG_H3K27ac_peaks_hmap_me3_Peaks.pdf --alpha 0.9 \
				--heatmapHeight 12 --heatmapWidth 6 \
					--refPointLabel "" -max 23

# H3K27me3 signal
computeMatrix reference-point -R H3K27me3.peaks -S BW_FILES/sorted_OLG_H3K27ac*bw -o OLG_H3K27ac_signal_me3_Peaks.matrix.gz -b 10000 -a 10000 -p 32
plotHeatmap -m OLG_H3K27ac_signal_me3_Peaks.matrix.gz -o OLG_H3K27ac_peaks_hmap_me3_Peaks.pdf --alpha 0.9 \
				--heatmapHeight 12 --heatmapWidth 6 \
					--refPointLabel "" -max 23



## compute signal for modality at peaks from BOTH modalities

## H3K27ac
computeMatrix reference-point \
	-R H3K27ac.peaks H3K27me3.peaks \
		-S BW_FILES/sorted_OLG_H3K27ac*bw \
			-o OLG_H3K27ac_signal.matrix.gz \
				--samplesLabel H3K27ac \
					-b 10000 -a 10000 -p 32

plotHeatmap -m OLG_H3K27ac_signal.matrix.gz \
	-o OLG_H3K27ac_peaks_hmap.pdf \
			--alpha 0.9 \
				--heatmapHeight 12 --heatmapWidth 6 \
					--refPointLabel "" \
						--samplesLabel "H3K27ac"


## H3K27me3
computeMatrix reference-point \
	-R H3K27me3.peaks H3K27ac.peaks \
		-S BW_FILES/sorted_OLG_H3K27me3*bw \
			-o OLG_H3K27me3_signal.matrix.gz \
				--samplesLabel H3K27me3 \
					-b 10000 -a 10000 -p 32

plotHeatmap -m OLG_H3K27me3_signal.matrix.gz \
	-o OLG_H3K27me3_peaks_hmap.pdf \
			--alpha 0.9 \
				--heatmapHeight 12 --heatmapWidth 6 \
					--refPointLabel "" \
						--samplesLabel "H3K27me3"