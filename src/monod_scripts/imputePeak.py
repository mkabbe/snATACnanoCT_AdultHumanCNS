import os
import logging
import anndata as ad
import scanpy as sc
import pandas as pd
import magic 
import scprep

logging.basicConfig(filename='logs/MAGIC_peakmatrix.log', format='%(asctime)s: %(message)s', level=logging.DEBUG)

#adata = ad.read("data/feature_matrices/210510_promoter_activity_matrix.h5ad")
adata = ad.read("data/feature_matrices/210510_LSI-2_RESULTS_merged_peak_matrix.h5ad")
#logging.info("Loaded GA matrix. Now creating MAGIC operator.")
logging.info("Loaded merged peak matrix. Now creating MAGIC operator.")

magic_op = magic.MAGIC(t=7)
logging.info("Created MAGIC operator. Now running MAGIC.")


adata_magic = magic_op.fit_transform(adata, genes="all_genes")
logging.info("Finished running MAGIC!")

logging.info("Converting MAGIC matrix to dataframe")
magic_df = pd.DataFrame(adata_magic.X, columns = adata_magic.var_names.tolist())

logging.info("Saving dataframe")
magic_df.to_csv("data/magic/210511_LSI-2_results_magic_matrix.csv", header=True, index_labels=False)
