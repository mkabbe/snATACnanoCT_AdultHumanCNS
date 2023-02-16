import os
import sys
import scanpy as sc
import anndata as ad
import numpy as np
import pandas as pd
import episcanpy.api as epi
import pysam
from scipy.sparse import lil_matrix
import argparse
from tqdm import tqdm
import pickle


## episcanpy environment should be activated

parser = argparse.ArgumentParser(description="Build custom feature matrix")
parser.add_argument("-f", "--fragments", help="path to fragments.tsv.gz file")
parser.add_argument("-p", "--peaks", help="path to the peak/features BED file")
parser.add_argument("-c", "--csv", help="path to the BARCODES PICKLE")
parser.add_argument("-o", "--out", help="name of the output .h5ad file containing the matrix")
parser.add_argument("-m", "--meta", help="metadata file")
parser.add_argument("-l", "--loaded", help="name of loaded output .h5ad")
args = parser.parse_args()

## check conda environment
#if sys.executable.split('/')[-3] == "episcanpy":
#    pass
#else:
#    print("need to activate conda environment\n")
#    sys.exit(1)

#check arguments
if not os.path.exists(args.fragments):
    print("fragments file not in the specified path\n")
    sys.exit(1)
if not os.path.exists(args.peaks):
    print("peaks file not in the specified path\n")
    sys.exit(1)
if not os.path.exists(args.csv):
    print("singlecell.csv file not in the specified path\n")
    sys.exit(1)


def build_mtx(tsv_file, annotation, barcodes, save=False):
    print("loading barcodes")
    
    # barcodes
    nb_barcodes = len(barcodes)
    dict_barcodes = {barcodes[i]: i for i in range(0,nb_barcodes)}
    set_barcode = set(barcodes)
    
    #load tabix
    tbx = pysam.TabixFile(tsv_file)
    
    window_list = []
    for chrom in sorted(annotation.keys()):
        window_list += [["".join(['chr', chrom]), int(n[0]), int(n[1])] for n in annotation[chrom]]
        
    print("building count matrix")
    mtx = lil_matrix((nb_barcodes, len(window_list)), dtype=np.float32)
    for i, tmp_feat in enumerate(tqdm(window_list)):
        for row in tbx.fetch(tmp_feat[0], tmp_feat[1], tmp_feat[2], parser=pysam.asTuple()):
            bcd = str(row).split('\t')[-2]
            if bcd in set_barcode:
                mtx[dict_barcodes[bcd], i] += 1
            else:
                continue
        
    print('building AnnData object')
    mtx = ad.AnnData(mtx.tocsr(), 
    obs=pd.DataFrame(index=barcodes), 
    var=pd.DataFrame(index=['_'.join([str(p) for p in n]) for n in window_list]))
    
    if save:
        mtx.write(save)
    
    return mtx

# load features
annot = epi.ct.load_features(args.peaks)

#load barcodes
bcd_ = pickle.load(open(args.csv, "rb"))

adata = build_mtx(tsv_file = args.fragments,
                    barcodes = bcd_, 
                    annotation = annot, 
                    save = args.out)


BASE = "/data/proj/GCB_MK/CZI/episcanpy_analysis/AGG_ATAC_210218/"
PICKLE = BASE + "ref/_sampleID.pickle"

epi.pp.filter_cells(adata, min_features=1)
epi.pp.filter_features(adata, min_cells=1)
adata.obs["log_nb_features"] = [np.log10(x) for x in adata.obs["nb_features"]]

with open(PICKLE, "rb") as f:
    GEM_TO_NGI_ID = pickle.load(f)
NGI_ID_TO_GEM = {v: k for k, v in GEM_TO_NGI_ID.items()}

bcd_list = adata.obs_names.to_list()
NGI_ID_list = [GEM_TO_NGI_ID[bcd.split("-")[-1]] for bcd in bcd_list]
adata.obs["NGI_ID"] = NGI_ID_list

df = pd.read_csv(args.meta, sep=",")
meta = df[df["NGI_ID"].notna()]

sample_attributes = meta.columns.to_list()
sample_attributes.remove("NGI_ID") #already added, not needed again

#populate anndata.obs
for attribute in sample_attributes:
    adata.obs[attribute] = [meta.loc[meta.NGI_ID == i, attribute].iloc[0] for i in NGI_ID_list]

#adata.write(BASE+"data/feature_matrices/210427_LOADED_2kb_matrix.h5ad")
adata.write(args.loaded)
