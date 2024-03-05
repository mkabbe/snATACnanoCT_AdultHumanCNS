import os
import logging
import anndata as ad
import scanpy as sc
import pandas as pd

logging.basicConfig(filename='logs/MAGIC_promotermatrix.log', format='%(asctime)s: %(message)s', level=logging.DEBUG)

#adata = ad.read("data/feature_matrices/210609_5kbpromoter_metagene_activity_matrix.h5ad")
#logging.info("Loaded GA matrix. Now creating MAGIC operator.")

#logging.info("Getting valid gene list for MAGIC")
#rna = ad.read("data/feature_matrices/HCA_RNA_all_annot.h5ad")
#rna.var_names = rna.var["gene_names"]
#var_names = rna.var_names.intersection(adata.var_names).tolist()

#logging.info(f"Running MAGIC now on {len(var_names)} genes")
#adata_magic = sc.external.pp.magic(adata, name_list=var_names, solver="approximate", knn_dist = "cosine", copy=True)
#logging.info("Finished running MAGIC. Writing results to file.")

#adata_magic.write("data/magic/210615_5kbpromoter_magic_matrix.h5ad")




#### 210621 ####

logging.info("Filtering out non OLG cells")
adata = ad.read("data/feature_matrices/210617_OLG_LSI-2_results_matrix.h5ad")
olg_bcd = adata.obs["barcode"].tolist()
del adata
ol_adata = ad.read("data/feature_matrices/210609_5kbpromoter_metagene_activity_matrix.h5ad")
ol_adata = ol_adata[olg_bcd,:]
logging.info("Loading RNA object to get valid genes")
adata_ref = ad.read("data/feature_matrices/HCA_oligo_and_opc_LS.h5ad")
adata_ref.var_names = adata_ref.var["gene_names"].tolist()
var_names = adata_ref.var_names.intersection(ol_adata.var_names).tolist()
logging.info(f"Identified {len(var_names)} genes")
del adata_ref
ol_adata.var = ol_adata.var.reset_index()
ol_adata.var_names = ol_adata.var.gene_name
ol_adata = ol_adata[:, var_names]
logging.info(f"data has {ol_adata.X.shape[0]} cells and {ol_adata.X.shape[1]} genes")

logging.info(f"Running MAGIC now on {len(var_names)} genes")
ol_magic = sc.external.pp.magic(ol_adata, name_list=var_names, knn_dist = "cosine", copy=True)
logging.info("Finished running MAGIC. Writing results to file.")

ol_magic.var_names = ol_magic.var.gene_name.tolist()
del ol_magic.var["gene_name"]
del ol_magic.var["level_0"]
del ol_magic.var["index"]

ol_magic.write("data/magic/210621_5kbpromoter_OLG_magic_matrix.h5ad")

