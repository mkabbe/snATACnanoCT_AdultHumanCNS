import os
import numpy as np
from motif_utils import *
import pickle

os.chdir("/data/proj/GCB_MK/CZI/episcanpy_analysis/AGG_ATAC_210218/data/motif/")

PEAK_FASTA = "210510_merged_louvain_cluster_peaks.fa"
PEAK_BED = "../../data/peaks/210510_merged_louvain_cluster_peaks.bed"
JASPAR_MOTIF_BASE_ORDER = ['A', 'C', 'G', 'T']

## Build list of GC distribution across all peaks
GCdist = []
peak_fa_list=[]
peak_name_list=[]
for header,seq in iter_fasta(PEAK_FASTA):
    peak_fa_list.append(seq) # store the FASTA sequences in a list
    peak_name_list.append(header) # store the FASTA headers in a list (chr:start-stop)
    GCdist.append(get_peak_GC_counts(seq, counts=False)) # get GC distribution


GCbounds = []
nbins = 25
LOW_GC = 0.33
HIGH_GC = 0.7

for n, gc in enumerate(np.percentile(GCdist, np.linspace(0, 100, nbins + 1, endpoint=True), interpolation='lower')):
    if n == 0 or n == nbins:
        GCbounds += [gc]
        continue
    if gc >= LOW_GC and gc < HIGH_GC:
        GCbounds += [gc]
GCbins = sorted(list(set(zip(GCbounds, GCbounds[1:]))))



GCdict = {}
base_counter = {}
for gc in GCbins:
    GCdict[gc] = {}
    GCdict[gc]['counter'] = Counter({base: 0 for base in JASPAR_MOTIF_BASE_ORDER})
    GCdict[gc]['peaks'] = []
    
for num, peak in enumerate(peak_fa_list):
    peakGC, base_counter = get_peak_GC_counts(peak)
    gc = findbin(peakGC, GCbins)
    GCdict[gc]['peaks'].append(num)
    GCdict[gc]['counter'] += Counter(base_counter)

for gc in GCbins:
        GCdict[gc]['total_bases'] = sum(GCdict[gc]['counter'].values())
        freq = {}
        for base in JASPAR_MOTIF_BASE_ORDER:
            freq[base] = np.divide(GCdict[gc]['counter'][base] + pseudocount,
                                   GCdict[gc]['total_bases'] + 4 * pseudocount, dtype=float)

        bg = {}
        bg['A'] = (freq['A'] + freq['T']) / 2
        bg['C'] = (freq['G'] + freq['C']) / 2
        bg['G'] = bg['C']
        bg['T'] = bg['A']
        GCdict[gc]['counter'] = [bg[base] for base in JASPAR_MOTIF_BASE_ORDER]


## write GC_bins information to file
GCdict_paths = {}
if not os.path.exists("GCdict_dump/"):
    os.mkdir("GCdict_dump/")
    
for gc in GCbins:
    GCdict_paths[gc] = f"GCdict_dump/GCdict_{gc[0]}_{gc[1]}"
    with open(GCdict_paths[gc], 'wb') as dump:
        pickle.dump(GCdict[gc], dump)

## MOODS needs to be run on each bin separately. 
chunk_bins = [{'skip': False, 'GCdict': GCdict_paths[chunk]} for chunk in GCbins]



    
