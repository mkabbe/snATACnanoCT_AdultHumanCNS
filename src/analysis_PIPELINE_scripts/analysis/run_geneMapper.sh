
#conda activate rose_py2

## required before running CRC analysis

#for CTYPE in AST MIGL CXEX CXINH MOL

for CTYPE in CXINH_ATAC #AST_ATAC #OPC_ATAC

do

	## connect enhancer to gene
#	python ROSE_geneMapper.py -i ${CTYPE}_roseSE/sorted_${CTYPE}_H3K27ac_peaks_AllEnhancers.table.txt -g HG38
	python ROSE_geneMapper.py -i ${CTYPE}_roseSE/sorted_${CTYPE}_fragments_AllEnhancers.table.txt -g HG38

	## replace header row of geneMapper output file with geneMap.header
	cp geneMap.header ${CTYPE}_ENHANCER_TO_GENE.txt
#	tail -n +2 ${CTYPE}_roseSE/sorted_${CTYPE}_H3K27ac_peaks_AllEnhancers_ENHANCER_TO_GENE.txt >> ${CTYPE}_ENHANCER_TO_GENE.txt
	tail -n +2 ${CTYPE}_roseSE/sorted_${CTYPE}_fragments_AllEnhancers_ENHANCER_TO_GENE.txt >> ${CTYPE}_ENHANCER_TO_GENE.txt
	mv ${CTYPE}_ENHANCER_TO_GENE.txt ${CTYPE}_roseSE/${CTYPE}_ENHANCER_TO_GENE.txt
done


