import os
from tqdm import tqdm
import anndata as ad


sample_list = []
for fname in os.listdir("../cellranger_COUNT_output"):
    sample_list.append(f"{fname}/outs/2kb_archrMatrix.h5ad")

ad_dict = {}
for count, sample in enumerate(sample_list):
    ad_dict[count] = ad.read(sample)

ad_list = []
for key in ad_dict:
    if key == 0:
        a = ad_dict[0]
    else:
        ad_list.append(ad_dict[key])

adata = a.concatenate(ad_list)
adata.write("/data/proj/GCB_MK/CZI/episcanpy_analysis/AGG_ATAC_210218/PIPELINE/data/matrices/211206_2kb_archr_matrix.h5ad")

