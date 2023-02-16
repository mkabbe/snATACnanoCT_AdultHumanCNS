import os 
import MOODS.parsers
import MOODS.tools
import MOODS.scan
import pyfaidx
from itertools import groupby
from collections import Counter
import random, string

def get_peak_GC_counts(peak, counts=True):
    '''Get GC% in base seq in a peak and (optionally) nucleotide counts'''

    seq = peak
    base_counter = {}
    base_counter['A'] = seq.count('A') + seq.count('a')
    base_counter['G'] = seq.count('G') + seq.count('g')
    base_counter['C'] = seq.count('C') + seq.count('c')
    base_counter['T'] = seq.count('T') + seq.count('t')
    base_counter['N'] = seq.count('N') + seq.count('n')
    peakGC = (base_counter['G'] + base_counter['C']) / sum(base_counter.values()) if len(seq) > 0 else 0.
    if counts:
        return peakGC, base_counter
    else:
        return peakGC
        
        

def findbin(peakGC, GCbins):
    ''''
    Identifies the bin each peak belongs to
    '''
    for gc in GCbins:
        if peakGC < gc[1] and peakGC >= gc[0]:
            return gc
    if abs(peakGC - GCbins[-1][1]) < 1e-5:
        return GCbins[-1]
        
        

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
    
    
    
def peak_reader(peak_fasta, select=[], only_list=False):
    '''
    First collects the peak headers and FASTA sequences in separate lists.
    If only_list = True, returns these 2 lists
    Else, generates the "bin-aware" peaks, which contains the peak list NOT in that GC bin.
    Returns the peaks and sequences as a pyfasta object.
    '''
    peak_fa_list=[]
    peak_name_list=[]

    for header,seq in iter_fasta(peak_fasta):
        peak_fa_list.append(seq) # store the FASTA sequences in a list
        peak_name_list.append(header) # store the FASTA headers in a list (chr:start-stop)

    if only_list:
        return(peak_fa_list, peak_name_list)

    peak_idx = set(select) # set of peaks to exclude
    peak_iter = []
    tmp_file = "TMP_FILES/" + ''.join(random.choice(string.ascii_uppercase + string.digits) for _ in range(6)) + ".tmp"
    for i in range(len(peak_fa_list)):
        if i in peak_idx:
            continue

        chrom, coords = peak_name_list[i].split(":")
        start, end = coords.split("-")

        peak_iter.append([chrom, int(start), int(end), "."]) # bed coords for the peak
        with open(tmp_file, "a") as outfa:
            outfa.write(f">{peak_name_list[i]}\n")
            outfa.write(f"{peak_fa_list[i]}\n")

    pyfaidx_obj = pyfaidx.Fasta(tmp_file)
    #os.system(f"rm {tmp_file}") # delete the tmp file

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
