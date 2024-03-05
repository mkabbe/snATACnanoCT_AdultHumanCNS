

DATADIR="/data/proj/GCB_MK/CZI/episcanpy_analysis/AGG_ATAC_210218/data/cluster_fragments/CELLTYPE_221210"

cd /data/proj/GCB_MK/scCT/nanoCT_EAE/rgt_analysis/ATAC_All

#mkdir -p Footprints
#mkdir -p MotifMatching
#mkdir -p DiffFootprinting


##	FOOTPRINTING
#for celltype in AST CBEX CBGRA CBINH CXEX CXINH ENDO MIGL MOL OPC;
#do
#	rgt-hint footprinting --atac-seq --paired-end --organism=hg38 \
#	--output-location=Footprints \
#	--output-prefix=${celltype} \
#	${DATADIR}/BAM_FILES/${celltype}.bed.bam \
#	${DATADIR}/PEAKS/${celltype}.bed_peaks.narrowPeak
#done



##	MOTIF MATCHING
#rgt-motifanalysis matching --organism=hg38 \
#--output-location=./MotifMatching \
#--input-files ./Footprints/*bed



##	DIFF FOOTPRINTING
rgt-hint differential --organism=hg38 --bc --nc 56 --window-size 500 \
--mpbs-files=./MotifMatching/AST_mpbs.bed,\
./MotifMatching/CBEX_mpbs.bed,\
./MotifMatching/CBGRA_mpbs.bed,\
./MotifMatching/CBINH_mpbs.bed,\
./MotifMatching/CXEX_mpbs.bed,\
./MotifMatching/CXINH_mpbs.bed,\
./MotifMatching/ENDO_mpbs.bed,\
./MotifMatching/MIGL_mpbs.bed,\
./MotifMatching/MOL_mpbs.bed,\
./MotifMatching/OPC_mpbs.bed \
--reads-files=${DATADIR}/BAM_FILES/AST.bed.bam,\
${DATADIR}/BAM_FILES/CBEX.bed.bam,\
${DATADIR}/BAM_FILES/CBGRA.bed.bam,\
${DATADIR}/BAM_FILES/CBINH.bed.bam,\
${DATADIR}/BAM_FILES/CXEX.bed.bam,\
${DATADIR}/BAM_FILES/CXINH.bed.bam,\
${DATADIR}/BAM_FILES/ENDO.bed.bam,\
${DATADIR}/BAM_FILES/MIGL.bed.bam,\
${DATADIR}/BAM_FILES/MOL.bed.bam,\
${DATADIR}/BAM_FILES/OPC.bed.bam \
--conditions=AST,CBEX,CBGRA,CBINH,CXEX,CXINH,ENDO,MIGL,MOL,OPC \
--output-location=./DiffFootprinting_new \
--output-prefix=ATAC_All_new


