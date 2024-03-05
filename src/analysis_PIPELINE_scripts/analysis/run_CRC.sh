
## conda activate episcanpy

CTYPE=$1

cd /data/proj/GCB_MK/scCT/rose_SE/rose/${CTYPE}_roseSE


crc -e ${CTYPE}_ENHANCER_TO_GENE.txt -g HG38 -c ../../../../CZI/episcanpy_analysis/AGG_ATAC_210218/ref/fasta/split -o crc_out_3 -n ${CTYPE}_CRC \
	-s /date/gcb/GCB_MK/ATAC_ctype_fragments_231219/NARROWPEAKS/sorted_${CTYPE}_fragments.bed_peaks.narrowPeak

## H3K27ac peaks
#-s ../../../nanoCT_EAE/fragments/P29054/ctype/H3K27ac/NARROWPEAKS/sorted_${CTYPE}_H3K27ac_peaks.bed

## ATAC peaks
#-s /date/gcb/GCB_MK/ATAC_ctype_fragments_231219/NARROWPEAKS/sorted_${CTYPE}_ATAC_fragments.bed_peaks.narrowPeak
