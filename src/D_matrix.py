#!/bin/python

import os 
import numpy as np
from scipy.sparse import lil_matrix


## first build a peak_dict and gene_dict
peaks = "merged_louvain_cluster_peaks.bed"  ## pre-sorted
genes = "gene_promoters.bed"    ## 2kb promoter centred on TSS
peak_dict = {}
gene_dict = {}

with open(peaks) as p:
    for line in p:
        line = line.strip().split("\t")
        peak_dict["_".join(line)] = line

with open(genes) as g:
    for line in g:
        line = line.strip().split("\t")
        gene_dict[line[-1]] = line[:-1]
        
n_peaks  = len(peak_dict)
n_genes = len(gene_dict)

## build zero matrix: rows = genes, col = peaks
D = lil_matrix((n_genes,n_peaks))

gene_index = -1
for gene in gene_dict.keys():
    gene_index += 1
    gene_chrom = gene_dict[gene][0]
    gene_start = int(gene_dict[gene][1])
    gene_end = int(gene_dict[gene][2])
    peak_index = -1
    for peak in peak_dict.keys():
        peak_index +=1
        peak_chrom = peak_dict[peak][0] 
        if peak_chrom != gene_chrom: ## no overlap if on different chromosomes 
            break
        else:
            peak_start = int(peak_dict[peak][1])
            peak_end = int(peak_dict[peak][2])
        if (peak_end <= gene_start): # peak is before the gene. Check next peak.
            continue
        elif (gene_end <= peak_start): # peak is after the gene. Check next gene.
            break
        else: #peak overlaps gene
            D[gene_index, peak_index] = 1

