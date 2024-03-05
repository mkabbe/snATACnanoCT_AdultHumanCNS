


DATADIR="/data/proj/GCB_MK/scCT/nanoCT_EAE/fragments/P29054/ctype/H3K27ac"

cd /data/proj/GCB_MK/scCT/nanoCT_EAE/rgt_analysis/H3K27ac

mkdir -p Footprints
mkdir -p MotifMatching
mkdir -p DiffFootprinting


##	FOOTPRINTING
for celltype in MIGL AST OPC MOL CXEX CXINH;
do
	rgt-hint footprinting --atac-seq --paired-end --organism=hg38 \
	--output-location=Footprints \
	--output-prefix=${celltype} \
	${DATADIR}/BW_FILES/${celltype}_H3K27ac_fragments.bed.sorted.bam \
	${DATADIR}/NARROWPEAKS/sorted_${celltype}_H3K27ac_fragments.bed_peaks.narrowPeak
done



##	MOTIF MATCHING
rgt-motifanalysis matching --organism=hg38 \
--output-location=./MotifMatching \
--input-files ./Footprints/*bed



##	DIFF FOOTPRINTING
rgt-hint differential --organism=hg38 --bc --nc 56 --window-size 500 \
--mpbs-files=./MotifMatching/MIGL_mpbs.bed,\
./MotifMatching/AST_mpbs.bed,\
./MotifMatching/OPC_mpbs.bed,\
./MotifMatching/MOL_mpbs.bed,\
./MotifMatching/CXEX_mpbs.bed,\
./MotifMatching/CXINH_mpbs.bed \
--reads-files=${DATADIR}/BW_FILES/MIGL_H3K27ac_fragments.bed.sorted.bam,\
${DATADIR}/BW_FILES/AST_H3K27ac_fragments.bed.sorted.bam,\
${DATADIR}/BW_FILES/OPC_H3K27ac_fragments.bed.sorted.bam,\
${DATADIR}/BW_FILES/MOL_H3K27ac_fragments.bed.sorted.bam,\
${DATADIR}/BW_FILES/CXEX_H3K27ac_fragments.bed.sorted.bam,\
${DATADIR}/BW_FILES/CXINH_H3K27ac_fragments.bed.sorted.bam \
--conditions=MIGL,AST,OPC,MOL,CXEX,CXINH \
--output-location=./DiffFootprinting_new \
--output-prefix=H3K27ac_All_new


