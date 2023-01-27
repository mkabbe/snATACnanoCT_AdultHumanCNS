
Steps for separating out modalities in nanoCT:

1) demultiplex modalities (multi-nanoCT)  -> run debarcode.py 
2) run cellranger-atac count

Steps for custom_cell_calling for scCT or nanoCT data: 

3) Store the possorted.bam in directory BAM/
4) Run call_peaks.sh with arguments
5) Run process_barcodes.sh with arguments
6) Run cell_calling.py or run the jupyter notebook (preferable)




