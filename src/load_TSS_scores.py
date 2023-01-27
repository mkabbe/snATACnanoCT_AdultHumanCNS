import os
import sys
import pickle
import pandas as pd
import anndata as ad



TSS_SCORES_DIR = sys.argv[1]
PICKLE = "ref/_sampleID.pickle"

tss_files = [os.path.join(TSS_SCORES_DIR,f) for f in os.listdir(TSS_SCORES_DIR) if f.endswith("_TSS_obs.csv")]

tss_dict = {}

with open(PICKLE, "rb") as f:
    GEM_TO_NGI_ID = pickle.load(f)
NGI_ID_TO_GEM = {v: k for k, v in GEM_TO_NGI_ID.items()}


for tssf in tss_files:
    ngi_id = tssf.split("/")[-1].split("_TSS")[0]
    df = pd.read_csv(tssf)
    for bcd in df.barcode:
        score = round(df.loc[df.barcode==bcd, "TSS_score"].iloc[0],2)
        bcd = bcd.split("-")[0]+"-"+NGI_ID_TO_GEM[ngi_id]
        tss_dict[bcd] = score

with open("ref/_TSSscores.pickle","wb") as handle:
    pickle.dump(tss_dict, handle)

adata = ad.read("data/matrices/211112_3X1K_LSI2_loaded.h5ad")

adata.obs["TSSe_score"] = [tss_dict[x] for x in adata.obs.barcode]