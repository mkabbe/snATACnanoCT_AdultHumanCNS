#!/bin/python

import os
import sys
import scanpy as sc
import anndata as ad
import numpy as np
import pandas as pd
import episcanpy.api as epi
import argparse
import pickle 

## episcanpy environment should be activated

parser = argparse.ArgumentParser(description="Build custom feature matrix")
parser.add_argument("-f", "--fragments", help="path to fragments.tsv.gz file")
parser.add_argument("-p", "--peaks", help="path to the peak/features BED file")
parser.add_argument("-c", "--csv", help="path to the FORMATTED singlecell.csv file.")
parser.add_argument("-o", "--out", help="name of the output .h5ad file containing the matrix")
parser.add_argument("-m", "--meta", help="metadata file")
parser.add_argument("-l", "--loaded", help="name of loaded output .h5ad")
args = parser.parse_args()

## check conda environment
#if sys.executable.split('/')[-3] == "episcanpy":
#    pass
#else:
#    print("need to activate conda environment\n")
#    sys.exit(1)
    
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
adata = epi.ct.bld_mtx_fly(tsv_file = args.fragments,
                            csv_file = args.csv,
                            annotation = annot,
                            save = args.out)

BASE = "/data/proj/GCB_MK/CZI/episcanpy_analysis/AGG_ATAC_210218/"
PICKLE = BASE + "ref/_sampleID.pickle"

epi.pp.filter_cells(adata, min_features=1)
epi.pp.filter_features(adata, min_cells=1)
adata.obs["log_nb_features"] = [np.log10(x) for x in adata.obs["nb_features"]]

with open(PICKLE, "rb") as f:
    GEM_TO_NGI_ID = pickle.load(f)
NGI_ID_TO_GEM = {v: k for k, v in GEM_TO_NGI_ID.items()}

bcd_list = adata.obs_names.to_list()
NGI_ID_list = [GEM_TO_NGI_ID[bcd.split("-")[-1]] for bcd in bcd_list] 
adata.obs["NGI_ID"] = NGI_ID_list

df = pd.read_csv(args.meta, sep=",")
meta = df[df["NGI_ID"].notna()]

sample_attributes = meta.columns.to_list()
sample_attributes.remove("NGI_ID") #already added, not needed again

#populate anndata.obs
for attribute in sample_attributes:
    adata.obs[attribute] = [meta.loc[meta.NGI_ID == i, attribute].iloc[0] for i in NGI_ID_list]

#adata.write(BASE+"data/feature_matrices/210427_LOADED_2kb_matrix.h5ad")
adata.write(args.loaded)
