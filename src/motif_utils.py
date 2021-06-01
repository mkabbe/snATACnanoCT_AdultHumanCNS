import os
import numpy as np
from itertools import groupby
from collections import Counter, namedtuple

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

Peak = namedtuple('Peak', ['chrom', 'start', 'end', 'strand'])
def peak_reader(peaks, select=[]):
    curr = 0
    with open(peaks, 'r') as peak_file:
        for line, row in enumerate(peak_file):
            chrom,coords = row.strip("\n").split(":")
            start,end = coords.split("-")
            if curr < len(select) and line == select[curr]:
                curr += 1
                continue
            yield Peak(chrom, int(start), int(end), '.')

