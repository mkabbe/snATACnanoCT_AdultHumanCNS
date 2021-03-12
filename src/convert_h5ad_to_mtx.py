#!/bin/python
import os
import pandas as pd
import anndata as ad
import argparse
from scipy.io import mmwrite


parser = argparse.ArgumentParser(description="convert h5ad feature matrix to mtx")
parser.add_argument("-h", "--h5ad", action="store")
parser.add_argument("-o", "--outdir", action="store")
parser.add_argument("-m", "--meta", action="store", default="False")
args = parser.parse_args()

# Takes anndata object and extract feature matrix and metadata
# stores peaks.bed, barcodes.tsv and matrix.mtx in the outdir folder


#os.chdir("/data/proj/GCB_MK/CZI/SCALE_analysis")
destination = args.outdir

if not os.path.exists(destination):
    os.mkdir(destination)

## read non-binary matrix
adata = ad.read(args.h5ad)


##read binary matrix
#adata =ad.read("210301_merged_peak_matrix.h5ad")

#write peaks
pd.DataFrame(adata.var.index).to_csv(os.path.join(destination, "peaks.bed"),   sep = "\t", index=False, header=False)
#write barcodes
pd.DataFrame(adata.obs.index).to_csv(os.path.join(destination, "barcodes.tsv"), sep = "\t", index=False, header=False)
#write matrix
mmwrite(os.path.join(destination, "matrix.mtx"), adata.X)


if args.meta:
    #write metadata
    adata.obs.to_csv(os.path.join(destination, "metadata.tsv"), sep="\t", index=True, header=True)
    
