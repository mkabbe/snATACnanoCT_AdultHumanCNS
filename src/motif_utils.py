import os
import numpy as np
from itertools import groupby
from collections import Counter, namedtuple
import pyfasta

## Functions adapted from 10X cellranger-atac

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



class Motifs:

    def __init__(self, dict_path, motif_path, bg=None, peaks=[]):
        
        self.dict_path = dict_path
        self.bg = bg
        self.peaks = peaks
        
        self.all_motifs = []
        for fn in os.listdir(motif_path):
            if fn.endswith(".jaspar"):
                self.all_motifs.append[motif_path.rstrip("/")+"/"+fn]

        with open(motif_path, "r") as infile:
            for line in infile:
                if line[0] == ">":
                    motif_id = line.lstrip(">").split()[0]
                    self.all_motifs.append(motif_id)
        
        
        self.genome_seq = pyfasta.Fasta("210510_merged_louvain_cluster_peaks.fa")


        
    def scan_motif_from_bed(self, peaks_iter, out_file=None, out_format="binary_bed", use_genome_bg=True, 
    pseudocount=1.0, pvalue=1e-7, window=7, tf_genes=None):
        if use_genome_bg:
            self.bg = [0.25,0.25,0.25,0.25]
        
        motif = self.get_motif_of(tf_genes)
        (matrices, thresholds) = self.prepare_moods_settings(motif_path=self.motif_path, bg = self.bg, pseudocount, pvalue)
        
        (seq, peaks_iter) = self.peak_reader(select=self.peaks) ## pyfasta object containing bin-specific peaks
         
        
        for peak_idx, peak in enumerate(peaks_iter):
            bed_coord = {'chr': peak[0],
                         'start': int(peak[1]),
                         'stop': int(peak[2]),
                         'name': peak_idx,
                         'strand': peak[3]}
        

    
    def get_motif_of(self, tf_genes):
        if tf_genes is None:
            return self.all_motifs
        selected = []
        for motif in self.all_motifs:
            if motif in tf_genes:
                selected.append(motif)
        return selected


    def prepare_moods_settings(self, motif_path, bg, pseudocount, pvalue=1e-7):
        
        matrices = [MOODS.parsers.pfm_to_log_odds(f, bg, pseudocount) for f in motif_path]
        matrices = matrices + [MOODS.tools.reverse_complement(m) for m in matrices]
        thresholds = [MOODS.tools.threshold_from_p(m, bg, pvalue) for m in matrices]
        
        return(matrices, thresholds)
        
    
        
    @staticmethod
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
    
    
    @staticmethod
    def peak_reader(select=[], only_list=False):
        '''
        First collects the peak headers and FASTA sequences in separate lists. 
        If only_list = True, returns these 2 lists
        Else, generates the "bin-aware" peaks, which contains the peak list NOT in that GC bin. 
        Returns the peaks and sequences as a pyfasta object.
        '''
        peak_fa_list=[]
        peak_name_list=[]
        
        for header,seq in iter_fasta("210510_merged_louvain_cluster_peaks.fa"):
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
        
        pyfasta_obj = pyfasta.Fasta("tmp", key_fn=lambda x: x.split()[0])
        os.system("rm tmp") # delete the tmp file
        
        return(pyfasta_obj, peak_iter)
        
