import os
import logging
import time
import anndata as ad
import numpy as np
import pandas as pd
import episcanpy.api as epi
import pysam
from scipy.sparse import lil_matrix
from tqdm import tqdm
import argparse
import pickle

## episcanpy environment should be activated

parser = argparse.ArgumentParser(description="Build custom feature matrix")
parser.add_argument("-f", "--fragments", help="path to fragments.tsv.gz file")
parser.add_argument("-p", "--peaks", help="path to the peak/features BED file")
parser.add_argument("-c", "--csv", help="path to the FORMATTED singlecell.csv file.")
parser.add_argument("-o", "--out", help="name of the output .h5ad file containing the matrix")
parser.add_argument("-m", "--meta", help="metadata file")
parser.add_argument("-r", "--results", help="name of loaded output .h5ad")
parser.add_argument("-l", "--log", help="name of log file")
parser.add_argument("-L", "--lsi", action="store_true", help="run default LSI?")
args = parser.parse_args()

logging.basicConfig(filename = args.log, format='%(asctime)s: %(message)s', level=logging.DEBUG)

init_time = time.time()

# load features
annot = epi.ct.load_features(args.peaks)

logging.info("Checking paths to files...")
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

logging.info("Files are OK. Proceeding to build matrix.")


## if the fragments file is too large to read into Pandas, the job fails. 
try:
    # build the count matrix. Takes a lot of time.
    adata = epi.ct.bld_mtx_fly(tsv_file = args.fragments,
                                csv_file = args.csv,
                                annotation = annot,
                                save = args.out)
except MemoryError:
    logging.info("Fragments file is too large, chunking and reading now.")
    
    def bld_chunked_mtx(tsv_file, annotation, csv_file=None, genome=None, save=False):
        print('loading barcodes')
        chunksize = 10**8
        barcodes = set()
        for chunk in pd.read_csv(tsv_file, sep='\t', header=None, chunksize=chunksize):
            chunked_bcd = sorted(chunk.loc[:,3].unique().tolist())
            barcodes = barcodes.union(chunked_bcd)

         # barcodes
        nb_barcodes = len(barcodes)
        dict_barcodes = {barcodes[i]: i for i in range(0, len(barcodes))}

        # Load tabix
        tbx = pysam.TabixFile(tsv_file)

            # format annotations
        window_list = []

        if genome:
            for chrom in sorted(annotation.keys()):
                window_list += [["".join([genome, '_chr', chrom]), int(n[0]), int(n[1])] for n in annotation[chrom]]
        else:
            for chrom in sorted(annotation.keys()):
                window_list += [["".join(['chr', chrom]), int(n[0]), int(n[1])] for n in annotation[chrom]]

        print('building count matrix')
        mtx = lil_matrix((nb_barcodes, len(window_list)), dtype=np.float32)
        for i, tmp_feat in enumerate(tqdm(window_list)):
            for row in tbx.fetch(tmp_feat[0], tmp_feat[1], tmp_feat[2], parser=pysam.asTuple()):
                mtx[dict_barcodes[str(row).split('\t')[-2]], i] += 1

        print('building AnnData object')
        mtx = ad.AnnData(mtx.tocsr(),
                            obs=pd.DataFrame(index=barcodes),
                            var=pd.DataFrame(index=['_'.join([str(p) for p in n]) for n in window_list]))

        if csv_file:
            print('filtering barcodes')
            df = pd.read_csv(csv_file)
            if genome == 'mm10':
                df_filtered = df[(df.is_mm10_cell_barcode == 1)]
            elif genome == 'hg19':
                df_filtered = df[(df.is_hg19_cell_barcode == 1)]
            else:
                df_filtered = df[(df.is_cell_barcode == 1)]

            barcodes = set(df_filtered.barcode.tolist())
            mtx = mtx[[i in barcodes for i in mtx.obs.index]].copy()

        if save:
            mtx.write(save)

        return mtx

    # build the count matrix. Takes a lot of time.
    adata = bld_chunked_mtx(tsv_file = args.fragments,
                                csv_file = args.csv,
                                annotation = annot,
                                save = args.out)

logging.info("Finished building matrix")

BASE = "/data/proj/GCB_MK/CZI/episcanpy_analysis/AGG_ATAC_210218/PIPELINE"
PICKLE = BASE + "ref/_sampleID.pickle"
#PICKLE = "ref/_sampleID.pickle"

epi.pp.filter_cells(adata, min_features=1)
epi.pp.filter_features(adata, min_cells=1)
adata.obs["log_nb_features"] = [np.log10(x) for x in adata.obs["nb_features"]]
adata.obs["barcode"] = adata.obs_names.copy()
logging.info("Filtered out emtpy cells and features")

with open(PICKLE, "rb") as f:
    GEM_TO_NGI_ID = pickle.load(f)
NGI_ID_TO_GEM = {v: k for k, v in GEM_TO_NGI_ID.items()}

bcd_list = adata.obs_names.to_list()
NGI_ID_list = [GEM_TO_NGI_ID[bcd.split("-")[-1]] for bcd in bcd_list]
adata.obs["NGI_ID"] = NGI_ID_list

logging.info("Filtered out emtpy cells and features")

#df = pd.read_csv(args.meta, sep=",")
#meta = df[df["NGI_ID"].notna()]

meta = pd.read_csv(args.meta, sep=",")
sample_attributes = meta.columns.to_list()
sample_attributes.remove("NGI_ID")
#populate anndata.obs
for attribute in sample_attributes:
    adata.obs[attribute] = [meta.loc[meta.NGI_ID == i, attribute].iloc[0] for i in adata.obs["NGI_ID"].tolist()]

if args.lsi:
    logging.info("Loaded sample attributes. Running LSI...")
    from lsi_functions import runTFIDF,runLSI   
    runLSI(adata)

else:
    logging.info("Loaded sample attributes. LSI will not be run.")

adata.write(args.results)
logging.info(f"Object saved to {args.results}")
logging.info(f"*** Finished building matrix in {np.round((time.time() - init_time)/60,2)} minutes *** \n")