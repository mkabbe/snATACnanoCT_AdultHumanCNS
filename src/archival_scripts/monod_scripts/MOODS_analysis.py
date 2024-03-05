import os
import sys
import MOODS.parsers
import MOODS.tools
import MOODS.scan
import pickle
import scipy
import anndata as ad
import episcanpy.api as epi
from motif_utils import *
from scipy.sparse import lil_matrix

os.chdir("/data/proj/GCB_MK/CZI/episcanpy_analysis/AGG_ATAC_210218/data/motif/")

PEAK_FASTA = sys.argv[1]
PEAK_BED = sys.argv[2]
PM_FILE = sys.argv[3]

## get list of all motifs
all_motifs = []
for fn in os.listdir("pfm_files/"):
            if fn.endswith(".jaspar"):
                all_motifs.append("pfm_files/"+fn)

## convert JASPAR to PFM 
for fn in all_motifs:
    if not os.path.exists(fn.split(".jaspar")[0]+".pfm"):
        jaspar_to_pfm(fn)

## get list of PFM motifs
pfm_motifs = []
for fn in os.listdir("pfm_files/"):
            if fn.endswith(".pfm"):
                pfm_motifs.append("pfm_files/"+fn)


## run MOODS for each GCbin 
for filename in os.listdir("GCdict_dump/"):
    
    ## load dict containing peak and background info for the GCbin
    bin_GCdict = pickle.load(open("GCdict_dump/"+filename, "rb"))
    
    bg = bin_GCdict["counter"] ## GC-aware background 
    peaks_in_bin = bin_GCdict["peaks"] ## indices for peaks in the bin
    pseudocount = 1.0

    (matrices, thresholds) = prepare_moods_settings(pfm_motifs, bg, pseudocount) ##load for MOODS analysis
    
    (peak_fa, peaks_iter) = peak_reader(peak_fasta=PEAK_FASTA, select=peaks_in_bin) ## load peaks in pyfaidx and get NON-binned peaks
    
    window = 7 ## size of the scanning window (in bp)
    scanner = MOODS.scan.Scanner(window)
    scanner.set_motifs(matrices, bg, thresholds)

    out_file = "moods_bin_results/moods_out_" + filename.split("GCdict_")[-1] + ".txt"
    out = open(out_file, "w")
        
    for peak_idx, peak in peak_fa.items():
    
        results = scanner.scan(str(peak))
        parsed = parse_scan_results(results, pfm_motifs, peak_idx)
            
        out.writelines(['\t'.join(map(str, item)) + '\n' for item in parsed])
        
    out.close()

MERGE_RESULTS = '''cat moods_bin_results/* | sort -k 1,1 -k2,2 | awk 'OFS="\t" {print($1"_"$2"_"$3,$4)}' > sorted_merged_results.txt'''
os.system(MERGE_RESULTS)


## store peaks
#peak_file = "../peaks/210510_merged_louvain_cluster_peaks.bed"
peak_file = PEAK_BED
peak_dict = {}
with open(peak_file, "r") as p:
    peak_idx = -1
    for line in p:
        peak_idx+=1
        line = "_".join(line.strip().split("\t")[:3])
        peak_dict[line] = peak_idx
        
## store peak_motif matches
peak_motif_dict = {}
motif_list = []
with open("sorted_merged_results.txt", "r") as infile:
    for line in infile:
        peak,motif = line.strip().split("\t")
        if peak not in peak_motif_dict.keys():
            peak_motif_dict[peak] = [motif]
        else:
            peak_motif_dict[peak].append(motif)
        
        if motif not in set(motif_list):
            motif_list.append(motif)

## store motif names
motif_dict = {}
motif_idx = -1
for motif in motif_list:
    motif_idx+=1
    motif_dict[motif] = motif_idx
    

## build PEAK x MOTIF matrix
PM = lil_matrix((len(peak_dict), len(motif_dict)))

for p in peak_motif_dict.keys():
    for m in peak_motif_dict[p]:
        PM[peak_dict[p], motif_dict[m]] += 1 ## not binary yet!

PM = PM.tocsr(copy=True) ## convert to CSR

## create AnnData object
pm_adata = ad.AnnData(PM)
pm_adata.obs_names = list(peak_dict.keys())
pm_adata.var_names = list(motif_dict.keys())

pm_adata.layers["raw_counts"] = pm_adata.X.copy() ## store raw counts
epi.pp.binarize(pm_adata)
pm_adata.layers["binarized_counts"] = pm_adata.X.copy() ## store binary counts

## save matrix
#pm_adata.write("../feature_matrices/210604_peak_by_motif_matrix.h5ad")
pm_adata.write(PM_FILE)
