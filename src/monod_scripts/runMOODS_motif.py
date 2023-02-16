import os 
from motif_utils import *
import sys
os.chdir("/data/proj/GCB_MK/CZI/episcanpy_analysis/AGG_ATAC_210218/data/motif/")

## get list of PFM motifs
pfm_motifs = []
for fn in os.listdir("pfm_files/"):
            if fn.endswith(".pfm"):
                pfm_motifs.append("pfm_files/"+fn)

for fn in os.listdir("GCdict_dump/"):
    bin_GCdict = pickle.load(open("GCdict_dump/"+fn, "rb"))
    bg = bin_GCdict["counter"]
    pseudocount = 1.0

    (matrices, thresholds) = prepare_moods_settings(pfm_motifs, bg, pseudocount)

    (peak_fa, peaks_iter) = peak_reader(select=bin_GCdict["peaks"])

    window = 7
    scanner = MOODS.scan.Scanner(window)  # parameter is the window size
    scanner.set_motifs(matrices, bg, thresholds)

    out_file = "moods_bin_results/moods_out_" + fn.split("GCdict_")[-1] + ".txt"
    out = open(out_file, "w")
    
    for peak_idx, peak in peak_fa.items(): 
    
        results = scanner.scan(str(peak))
        parsed = parse_scan_results(results, pfm_motifs, peak_idx)
            
        out.writelines(['\t'.join(map(str, item)) + '\n' for item in parsed])
        
    out.close()

