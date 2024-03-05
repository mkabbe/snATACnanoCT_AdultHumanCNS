import os
import numpy as np
import pandas as pd
import anndata as ad
import scanpy as sc
import episcanpy.api as epi
import magic
import scipy

os.chdir("/data/proj/GCB_MK/CZI/episcanpy_analysis/AGG_ATAC_210218")

gadata = ad.read("data/feature_matrices/211102_metagene_promoter_activity_matrix_3x10K.h5ad")
rna = ad.read("data/feature_matrices/HCA_RNA_all_annot.h5ad")
rna.var_names = rna.var.gene_names.copy()
markers = list(rna.var_names.intersection(gadata.var_names))
#markers = ['MOG','PLP1','SOX10','PDGFRA','OLIG2','GAD1','AQP4','TUBB3','AIF1','GFAP','RBFOX3','SLC17A7']

magic_op = magic.MAGIC(t=3)
gene_adata_magic = magic_op.fit_transform(gadata, genes=markers)

try:
    gene_adata_magic.X =  scipy.sparse.csr_matrix(gene_adata_magic.X)
    gene_adata_magic.write("data/feature_matrices/211109_metagene_promoter_activity_MAGIC_matrix_3x10K.h5ad")
except:
    gene_adata_magic.write("data/feature_matrices/211109_metagene_promoter_activity_MAGIC_matrix_3x10K.h5ad")
