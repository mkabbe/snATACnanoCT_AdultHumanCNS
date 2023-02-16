#!/bin/python
import os
import logging
import time
import csv
import scipy
import scanpy as sc
import anndata as ad
import numpy as np
import pandas as pd
import episcanpy.api as epi
from metagene_func import *


os.chdir("/data/proj/GCB_MK/CZI/episcanpy_analysis/AGG_ATAC_210218")
logging.basicConfig(filename='logs/MAGIC2mGA.log', format='%(asctime)s: %(message)s', level=logging.DEBUG)

init_time = time.time()

## Load MAGIC matrix
adata = ad.read("data/magic/210527_gene_activity_magic_matrix.h5ad")
logging.info("Loaded MAGIC matrix. \n")


## Load RNA data
logging.info("Reading RNA matrix and extracting marker genes..... \n")
rna = ad.read("data/feature_matrices/HCA_RNA_all_annot.h5ad")
gene_modules  = extractMarkerGenes(rna_adata = rna)
logging.info("Marker genes extracted! \n")

## Get metagene scores
logging.info("Computing metagene scores..... \n")
addMetageneScores(adata=adata, gene_modules=gene_modules)
logging.info("***Finished!*** \n")

adata.write("data/magic/210527_gene_metagene_activity_magic_matrix.h5ad")

logging.info("Total time elapsed: {} minutes. \n".format((time.time() - init_time)/60))
