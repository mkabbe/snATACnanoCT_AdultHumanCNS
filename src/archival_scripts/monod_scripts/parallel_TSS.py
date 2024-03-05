import os
import sys
import numpy as np
import pybedtools
import metaseq
from joblib import Parallel,delayed #,parallel_backend
import multiprocessing 

## Adapted from ENCODE-atac-seq-pipeline

def scoreTSS(bam_file, tss="ref/GRCh38.p13_TSS.bed", prefix=10, chromsizes="ref/hg38.chrom.sizes", read_len=51, bins=400, bp_edge=2000, processes=56, greenleaf_norm=True):
    # Load the TSS file
    barcode = bam_file.split("/")[-1].split("_")[0]
    tss = pybedtools.BedTool(tss)
    tss_ext = tss.slop(b=bp_edge, g=chromsizes)

    # Load the bam file
    bam = metaseq.genomic_signal(bam_file, 'bam')
    bam_array = bam.array(tss_ext, bins=bins, shift_width = -read_len/2, # Shift to center the read on the cut site
                          processes=processes, stranded=True)

    if greenleaf_norm:
        # Use enough bins to cover 100 bp on either end
        num_edge_bins = int(100/(2*bp_edge/bins))
        bin_means = bam_array.mean(axis=0)
        avg_noise = (sum(bin_means[:num_edge_bins]) +
                     sum(bin_means[-num_edge_bins:]))/(2*num_edge_bins)
        bam_array /= avg_noise
    else:
        bam_array /= bam.mapped_read_count() / 1e6

    tss_point_val = max(bam_array.mean(axis=0))
    with open("data/tss_csv/{}_TSSscore.csv".format(barcode), "a") as f:
        f.write(barcode + "," + str(tss_point_val))
#    return (barcode + "," + str(tss_point_val))

#processes = multiprocessing.cpu_count()
directory="data/P20602_1005_bamFiles/"
file_list = [x for x in os.listdir(directory) if x.endswith(".bam")]

for bam in file_list:
    scoreTSS(directory+bam)


#Parallel(n_jobs=6)(delayed(scoreTSS(directory+bam)) for bam in file_list)

#Parallel(n_jobs=6)(delayed(scoreTSS(directory+bam)) for bam in os.listdir(directory) if bam.endswith(".bam"))

#with parallel_backend("threading", n_jobs=12):
#    Parallel()(delayed(scoreTSS(directory+bam)) for bam in os.listdir(directory) if bam.endswith(".bam"))

