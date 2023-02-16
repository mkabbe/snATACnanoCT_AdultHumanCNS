import os
import logging
import anndata as ad
import scanpy as sc
import pandas as pd
import magic 
import scipy
import scprep
from metagene_func import extractMarkerGenes, addMetageneScores

logging.basicConfig(filename='logs/MAGIC_promotermatrix.log', format='%(asctime)s: %(message)s', level=logging.DEBUG)

#adata = ad.read("data/feature_matrices/210521_gene_activity_matrix.h5ad")
adata = ad.read("data/feature_matrices/210609_5kbpromoter_metagene_activity_matrix.h5ad")
logging.info("Loaded GA matrix. Now creating MAGIC operator.")

magic_op = magic.MAGIC(t=7, n_jobs=-2, solver="approximate", n_pca=30)
logging.info("Created MAGIC operator. Now running MAGIC.")


adata_magic = magic_op.fit_transform(adata, genes="all_genes")
logging.info("Finished running MAGIC!")

logging.info("Saving as dataframe...")
magic_df = pd.DataFrame(adata_magic.X, columns = adata_magic.var_names.tolist())
magic_df.to_csv("data/magic/210611_5kbpromoter_activity_magic_matrix.csv", header=True)#, index_labels=False)


logging.info("Converting MAGIC matrix to scipy.sparse csr_matrix")
adata_magic.X = scipy.sparse.csr_matrix(adata_magic.X)
logging.info("Saving h5ad")
adata_magic.write("data/magic/210611_5kbpromoter_activity_magic_matrix.h5ad")

## Load RNA data
logging.info("Reading RNA matrix and extracting marker genes..... \n")
rna = ad.read("data/feature_matrices/HCA_RNA_all_annot.h5ad")
gene_modules  = extractMarkerGenes(rna_adata = rna)
logging.info("Marker genes extracted! \n")

## Get metagene scores
logging.info("Computing metagene scores..... \n")
addMetageneScores(adata=adata_magic, gene_modules=gene_modules)
logging.info("***Finished!*** \n")
adata_magic.write("data/magic/210611_5kbpromoter_metagene_activity_magic_matrix.h5ad")

#adata_magic.write("data/magic/210524_gene_activity_magic_matrix.h5ad")
#magic_df.to_csv("data/magic/210514_promoter_activity_magic_matrix.csv", header=True)#, index_labels=False)
