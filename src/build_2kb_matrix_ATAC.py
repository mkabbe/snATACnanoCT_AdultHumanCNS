#!/bin/python

import os
import csv
import scanpy as sc
import anndata as ad
import numpy as np
import pandas as pd
import episcanpy.api as epi
import scipy


os.chdir("/Volumes/shared/Molekylär Neurobiologi/Castelo-Branco/Mukund/CZI_ADULT/CZI_ATAC")

RESULTS_FILE = "data/epi/Binned_2kb_count_AGG_ATAC_results.h5ad"
DATA_PATH = "../../../NGSDATA/scATAC_human_adult_CZI_1/Processed_data/AGG_ATAC_210218/outs/"
BINS_PATH = "bin/2kb_binned_hg38.bed"
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


annot = epi.ct.load_features(BINS_PATH)

## build the matrix (takes a *LOT* of time, approx. 3 hours)
counts = epi.ct.bld_mtx_fly(tsv_file = FRAG_PATH,
                            csv_file = CSV_PATH,
                           annotation = annot,
                           save = "data/epi/binned_2kb_count_matrix.h5ad")