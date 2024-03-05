#!/bin/python

import os
import csv
import scanpy as sc
import anndata as ad
import numpy as np
import pandas as pd
import episcanpy.api as epi

os.chdir("/Volumes/shared/Molekylär Neurobiologi/Castelo-Branco/Mukund/CZI_ADULT/CZI_ATAC")

RESULTS_FILE = "data/epi/peak_count_AGG_ATAC_results.h5ad"
DATA_PATH = "../../../NGSDATA/scATAC_human_adult_CZI_1/Processed_data/AGG_ATAC_210218/outs/"
PEAK_PATH = DATA_PATH + "peaks.bed"
FRAG_PATH = DATA_PATH + "fragments.tsv.gz"
CSV_PATH = DATA_PATH + "singlecell_EPI.csv"
METADATA_PATH = "~/Documents/PhD/snATAC-seq/CZI/Metadata/sample_randomizer_metadata.csv"
GEM_TO_NGI_ID = {"1" : "P20056_1001",
                "2" : "P20056_1002",
                "3" : "P20056_1003",
                "4" : "P20056_1004", 
                "5" : "P20057_1001",
                "6" : "P20057_1002",
                "7" : "P20057_1003"}
NGI_ID_TO_GEM = {v: k for k, v in GEM_TO_NGI_ID.items()}


#annot = epi.ct.load_features(PEAK_PATH)

## build the matrix (takes time)
#counts = epi.ct.bld_mtx_fly(tsv_file = FRAG_PATH,
#                            csv_file = CSV_PATH,
#                           annotation = annot,
#                           save = "data/epi/peak_count_matrix.h5ad")

adata = ad.read("data/epi/peak_count_matrix.h5ad")

bcd_list = adata.obs_names.to_list()
NGI_ID_list = [GEM_TO_NGI_ID[bcd.split("-")[-1]] for bcd in bcd_list] 
adata.obs["NGI_ID"] = NGI_ID_list

df = pd.read_csv(METADATA_PATH, sep=",")
meta = df[df["NGI_ID"].notna()]

sample_attributes = meta.columns.to_list()
sample_attributes.remove("NGI_ID") #already added, not needed again

#populate anndata.obs
for attribute in sample_attributes:
    adata.obs[attribute] = [meta.loc[meta.NGI_ID == i, attribute].iloc[0] for i in NGI_ID_list]


## get barcodes of doublets identified from ArchR.
ARCHR_ALL_BCD = "./archr_all_bcd.csv"
ARCHR_SINGLET_BCD = "./archr_singlet_bcd.csv"

with open(ARCHR_ALL_BCD, newline='') as f:
    reader = csv.reader(f)
    archr_all = list(reader)

with open(ARCHR_SINGLET_BCD, newline='') as f:
    reader = csv.reader(f)
    archr_singlets = list(reader)

#change archr_barcode naming style to match anndata.obs
archr_list = [archr_all, archr_singlets]
idx = 0
for ls in archr_list:
    tmp = []
    for archr_bcd in ls[1:]:
        ngi_id,bcd = archr_bcd[-1].split("#")
        tmp.append(bcd[:-1] + NGI_ID_TO_GEM[ngi_id])
    idx += 1
    if idx == 1:
        archr_all = tmp
    else:
        archr_singlets = tmp

# get doublet barcodes
archr_doublets = [bcd for bcd in archr_all if bcd not in archr_singlets]
#add doublet attribute
adata.obs["ArchR_doublet"] = ["True" if bcd in archr_doublets  else "False" for bcd in bcd_list]


#get gene annotations
epi.tl.find_genes(adata,
           gtf_file='ref/gencode.v35.annotation.gtf',
           key_added='transcript_annotation',
           upstream=2000,
           feature_type='transcript',
           annotation='HAVANA',
           raw=False)
           

# get "commonness" of each feature (make sure matrix is binarized)
n_cells = np.sum(adata.X, axis=0).tolist()
adata.var["n_cells"] = n_cells[0]

#save
adata.raw = adata
adata.write(RESULTS_FILE)


# For downstream analysis, load the RESULTS_FILE to fire up the populated anndata matrix

