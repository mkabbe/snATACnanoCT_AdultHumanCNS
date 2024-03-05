#!/bin/python

import os
import sys
import scanpy as sc
import anndata as ad
import numpy as np
import pandas as pd
import episcanpy.api as epi
import warnings
import argparse
from aux import *
warnings.filterwarnings("ignore")

## episcanpy environment should be activated
## parse arguments
parser = argparse.ArgumentParser(description="Build custom feature matrix")
parser.add_argument("-f", "--fragments", help="path to fragments.tsv.gz file")
parser.add_argument("-p", "--peaks", help="path to the peak/features BED file")
parser.add_argument("-c", "--csv", help="path to the FORMATTED singlecell.csv file.")
parser.add_argument("-o", "--out", help="name of the output .h5ad file containing the matrix")
parser.add_argument("-a", "--analyzed", help="name of the analyzed .h5ad file")
parser.add_argument("--clusters", help="barcode cluster file")
parser.add_argument("--clusterDir", help="directory for storing cluster fragments")
args = parser.parse_args()

## check for correct singlecell.csv file
if not os.path.isfile(os.getcwd()+ "/" + args.csv):
    os.system("sed 's/is__cell/is_cell/' AGG_ATAC_2104/outs/singlecell.csv > {}".format(args.csv))

## check conda environment
if sys.executable.split('/')[-3] == "episcanpy":
    pass
else:
    print("activate conda environment\n")
    sys.exit(1)
    
#check arguments
if not os.path.exists(args.fragments):
    print("fragments file not in the specified path\n")
    sys.exit(1)
if not os.path.exists(args.peaks):
    print("peaks file not in the specified path\n")
    sys.exit(1)
if not os.path.exists(args.csv):
    print("singlecell.csv file not in the specified path\n")
    sys.exit(1)
    
# load features
annot = epi.ct.load_features(args.peaks)

# build the count matrix. Takes a lot of time. 
counts = epi.ct.bld_mtx_fly(tsv_file = args.fragments,
                            csv_file = args.csv,
                            annotation = annot,
                            save = args.out)

## load anndata file
adata = ad.read(args.out)
if "barcode" not in adata.obs.columns():
    adata.obs["barcode"] = adata.obs_names

## pre-processing 
epi.pp.binarize(adata)
epi.pp.filter_cells(adata, min_features=1)
epi.pp.filter_features(adata, min_cells=1)
adata.obs["log_nb_features"] = [np.log10(x) for x in adata.obs["nb_features"]]
## filter out lowest 10% of cells
min_features = np.quantile(adata.obs["nb_features"], 0.1)
epi.pp.filter_cells(adata, min_features=min_features)
## select top 20,000 features
min_cells = np.sort(adata.var["n_cells"])[-20000]
adata = adata[:, adata.var["n_cells"] >= min_cells]

## TF-IDF normalization 
runTFIDF(adata)

## run SVD
sc.tl.pca(adata, zero_center=False, n_comps=51)

## remove first PC (depth correlated)
removefirstPC(adata)

## get neighbors, run UMAP and louvain clustering
sc.pp.neighbors(adata, n_pcs = 50, n_neighbors = 15, metric="cosine")
sc.tl.umap(adata, n_components=2)
epi.tl.louvain(adata)

## extract cluster_IDs 
bc_cluster_df = adata.obs[["barcode","louvain"]]
bc_cluster_df.to_csv(args.clusters, header=False, index=False)

## write to results file
adata.write(args.analyzed)

## split fragments by cluster
barcodes = parse_barcodes(args.clusters)
reads_dic = filter_fragments_by_cluster(args.fragments, barcodes)

for key in reads_dic:
    with open(args.clusterDir + "_" + key + ".bed", "w") as f:
        f.write("\n".join(reads_dic[key]))


### add MACS peak calling


## > build Agg template
## > run cellranger Agg
## > modify the singlecell.csv file
## > build the 2kb count matrix
## > LSI-1 analysis 
## > clustering + fragment filtering
## > call and merge peaks
## > build peak matrix
## > annotate cells and features


