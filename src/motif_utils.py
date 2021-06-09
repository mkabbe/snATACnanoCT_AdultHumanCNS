import os 
import MOODS.parsers
import MOODS.tools
import MOODS.scan
import pyfaidx
from itertools import groupby



def jaspar_to_pfm(filename):
        '''
        Convert the jaspar files to PFM format (MOODS compatible)
        ** Assumes filename ends in .jaspar
        '''
        outfile = filename.split(".jaspar")[0] + ".pfm"
        
        with open(filename, "r") as f:
            for line in f:
                if line.startswith(">"):
                    continue
                with open(outfile, "a") as out_:
                    out_.write("\t".join(line.lstrip("AGCT").strip().strip("\[\]").split())+"\n")



def iter_fasta(filename):
    ''' 
    returns generator object
    '''
    with open(filename, "r") as f:
        iterator = groupby(f, lambda line: line[0] == ">")
        for is_header, group in iterator:
            if is_header:
                header = group.__next__()[1:].strip()
            else: 
                yield header, "".join(s.strip() for s in group)



def prepare_moods_settings(motif_path, bg, pseudocount, pvalue=1e-7):
        
    matrices = [MOODS.parsers.pfm_to_log_odds(f, bg, pseudocount) for f in motif_path]
    matrices = matrices + [MOODS.tools.reverse_complement(m) for m in matrices]
    thresholds = [MOODS.tools.threshold_from_p(m, bg, pvalue) for m in matrices]
     
    return(matrices, thresholds)
    
    
    
def peak_reader(select=[], only_list=False):
    '''
    First collects the peak headers and FASTA sequences in separate lists. 
    If only_list = True, returns these 2 lists
    Else, generates the "bin-aware" peaks, which contains the peak list NOT in that GC bin. 
    Returns the peaks and sequences as a pyfasta object.
    '''
    peak_fa_list=[]
    peak_name_list=[]
    
    for header,seq in iter_fasta("../210510_merged_louvain_cluster_peaks.fa"):
        peak_fa_list.append(seq) # store the FASTA sequences in a list
        peak_name_list.append(header) # store the FASTA headers in a list (chr:start-stop)

    if only_list:
        return(peak_fa_list, peak_name_list)
    
    peak_idx = set(select) # set of peaks to exclude
    peak_iter = []
    for i in range(len(peak_fa_list)):
        if i in peak_idx:
            continue
        
        chrom, coords = peak_name_list[i].split(":")
        start, end = coords.split("-")
        
        peak_iter.append([chrom, int(start), int(end), "."]) # bed coords for the peak
        
        with open("tmp", "a") as outfa:
            outfa.write(f">{peak_name_list[i]}\n")
            outfa.write(f"{peak_fa_list[i]}\n")
        
    pyfaidx_obj = pyfaidx.Fasta("tmp")
    os.system("rm tmp") # delete the tmp file
     
    return(pyfaidx_obj, peak_iter)
    


def parse_scan_results(moods_scan_res, motifs, idx, out_format="binary-bed"):
    
    chrom,coords = idx.split(":")
    start,stop = coords.split("-")
    all_hits = []
    for (motif_idx, hits) in enumerate(moods_scan_res):
        motif = motifs[motif_idx % len(motifs)]
        strand = "-" if motif_idx >= len(motifs) else "+"
        motif_name = motif.split("/")[-1].split(".pfm")[0]
    
        if len(hits) > 0: 
            record = [chrom, int(start), int(stop), motif_name]
            all_hits.append(record)
            continue
    return all_hits