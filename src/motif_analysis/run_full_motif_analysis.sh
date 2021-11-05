#/bin/bash


# activate Environment
conda activate episcanpy

cd /data/proj/GCB_MK/CZI/episcanpy_analysis/AGG_ATAC_210218/data/motif

## store names (change accordingly)
PEAKS_BED= "../../data/peaks/211029_merged_leiden_cluster_peaks_3X10K.bed"
PEAKS_FASTA="211029_merged_leiden_cluster_peaks_3X10K.fa"
PM_MATRIX="../../data/feature_matrices/211105_3X10K_peak_by_motif_matrix.h5ad"

## not sure why but condor has issues recognizing python, so run this:
PYTHON="/home/mukund/miniconda3/envs/episcanpy/bin/python"


## download the human reference FASTA
if [[ ! -f ../../ref/fasta/Homo_sapiens.GRCh38.dna.primary_assembly.fa ]]; then
	
	wget http://ftp.ensembl.org/pub/release-104/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
	## unzip and move to ref directory	
	gunzip Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
	mv Homo_sapiens.GRCh38.dna.primary_assembly.fa ../../ref/fasta/

fi

## pre-process
sed -i 's/ dna.*$//g' ../../ref/fasta/Homo_sapiens.GRCh38.dna.primary_assembly.fa 
sed -i 's/>/>chr/g' ../../ref/fasta/Homo_sapiens.GRCh38.dna.primary_assembly.fa

## convert peaks BED to FASTA
bedtools getfasta -fi ../../ref/fasta/Homo_sapiens.GRCh38.dna.primary_assembly.fa -bed $PEAKS_BED -fo $PEAKS_FASTA

## prep GC bias bins 
$PYTHON ../../src/peak_GC_distribution.py $PEAKS_FASTA $PEAKS_BED

## Run MOODS and build Peak-Motif matrix
$PYTHON ../../src/MOODS_analysis.py $PEAKS_FASTA $PEAKS_BED $PM_MATRIX
