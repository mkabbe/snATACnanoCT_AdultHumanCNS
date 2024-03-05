import pandas as pd
import numpy as np
import anndata as ad
import scipy 
import logging 

logging.basicConfig(filename='logs/MAGIC_promotermatrix.log', format='%(asctime)s: %(message)s', level=logging.DEBUG)
adata = ad.read("data/feature_matrices/210510_LSI-2_RESULTS_merged_peak_matrix.h5ad")
logging.info("Loaded anndata object")


magic_df = pd.read_csv("data/magic/210514_promoter_activity_magic_matrix.csv")
logging.info("Loaded MAGIC dataframe")

magic_X = scipy.sparse.csr.csr_matrix(magic_df.values)
logging.info("Converted MAGIC to sparse matrix")

adata.X = magic_X
logging.info("Loaded MAGIC into anndata")

adata.write("data/magic/210524_promoter_magic_matrix.h5ad")
logging.info("Saved MAGIC anndata")
#try:
#	adata.X = magic_X
#	adata.write("data/magic/210524_promoter_magic_matrix.h5ad")
#except:
#	scipy.sparse.save_npz("data/magic/210524_promoter_magic_matrix.npz", magic_X)


