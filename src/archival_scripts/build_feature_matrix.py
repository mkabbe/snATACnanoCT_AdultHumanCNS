#!/bin/python
import os 
import scanpy as sc
import anndata as ad
import episcanpy.api as epi
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from bin_genome import bin_genome

# not sure if this part is only jupyter-compatible
plt.rcParam['figure.figsize']=[8,8]
%config InlineBackend.figure_format = 'retina'
!conda activate episcanpy

#generate features as windows spanning genome
binsize = 5000
bin_abbr = str(binsize//1000)+"kb"

if not os.path.isfile("ref/{}_bins_hg38.bed".format(bin_abbr):
    bin_genome(
        binsize=binsize,
        stepsize=0,
        chromsize_file="ref/hg38.chrom.sizes",
        outfile="ref/{}_bins_hg38.bed".format(bin_abbr)
    )
                   

#load features
features = epi.ct.load_features("ref/{}_bins_hg38.bed".format(bin_abbr))

#build feature matrix
if not os.path.isfile("data/epi/{}.feature_matrix.h5ad".format(bin_abbr)):
    counts = epi.ct.bld_mtx_fly(
        tsv_file = "data/raw/fragments.tsv.gz",
        csv_file = "data/raw/singlecell.csv",
        annotation = features,
        save = "data/epi/{}.feature_matrix.h5ad".format(bin_abbr)
    )

else:
    print("Feature matrix already exists.")
