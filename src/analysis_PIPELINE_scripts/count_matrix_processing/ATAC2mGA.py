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
logging.basicConfig(filename='logs/ATAC2mGA2.log', format='%(asctime)s: %(message)s', level=logging.DEBUG)

init_time = time.time()
 
## Load 2kb matrix
#adata = ad.read("data/feature_matrices/210921_3x1K_2kb_feature_matrix_LOADED.h5ad")
adata = ad.read("data/feature_matrices/210921_3x10K_2kb_feature_matrix.h5ad")
logging.info("Loaded ATAC matrix. \n")

## Extract the gene coordinates (gene body + 2kb upstream)
logging.info("Extracting gene promoter coordinates..... \n")
#logging.info("Extracting gene BODY + PROMOTER coordinates..... \n")
#gtf = extractGeneAnnot(gtf_file="ref/gencode.v35.annotation.gtf")
gtf = extractPromoterCoord(gtf_file="ref/gencode.v35.annotation.gtf", window = 5000)
#gtf = extractGeneAnnot(gtf_file="ref/gencode.v35.annotation.gtf", upstream = 2000)
logging.info("Extracted 5kb promoter coordinates! \n")
#logging.info("Extracted gene BODY + PROMOTER  coordinates! \n")
adata_features = extractFeatureCoordinates(adata)
logging.info("Extracted peak coordinates! \n")

## Build GA matrix
start = time.time()
logging.info("Building gene activity matrix....... \n")
gene_adata = buildGeneActivityMatrix(adata=adata, gtf=gtf, raw_adata_features=adata_features, feature_type="gene")
logging.info("Built gene activity matrix in {} minutes. \n".format((time.time() - start)/60))

logging.info("Converting to sparse matrix... \n")
gene_adata.X = scipy.sparse.csr_matrix(gene_adata.X)

logging.info("Saving... \n")

#gene_adata.write("data/feature_matrices/211101_gene_activity_from2kbmtx_all_cells.h5ad")
#gene_adata.write("data/feature_matrices/211102_gene_activity_from2kbmtx_3x10K.h5ad")
gene_adata.write("data/feature_matrices/211102_promoter_activity_from2kbmtx_3x10K.h5ad")

## Load RNA data
logging.info("Reading RNA matrix and extracting marker genes..... \n")
rna = ad.read("data/feature_matrices/HCA_RNA_all_annot.h5ad")
gene_modules  = extractMarkerGenes(rna_adata = rna)
logging.info("Marker genes extracted! \n")

## Get metagene scores
logging.info("Computing metagene scores..... \n")
addMetageneScores(adata=gene_adata, gene_modules=gene_modules)
logging.info("***Finished!*** \n")

#gene_adata.write("data/feature_matrices/211101_metagene_activity_matrix_all_cells.h5ad")
#gene_adata.write("data/feature_matrices/211102_metagene_activity_matrix_3x10K.h5ad")
gene_adata.write("data/feature_matrices/211102_metagene_promoter_activity_matrix_3x10K.h5ad")

logging.info("Total time elapsed: {} minutes. \n".format((time.time() - init_time)/60))

