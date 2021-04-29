import pandas as pd
import numpy as np
import pysam
from scipy.sparse import lil_matrix
import pickle

HUMAN = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', 
        '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', 
        '21', '22','X', 'Y']

def load_features(file_features, chromosomes=HUMAN, path=""):
    """
    load TSS features when input file is bed. 
    Returns dictionary with key=chromosome and values= list of features (each feature as list)
    """
    features_chrom = {}
    for c in chromosomes:
        features_chrom[c] = []
    with open(path + file_features) as features:
        for line in features:
            line = line.strip().split("\t")
            if line[0][3:] in chromosomes:
                try:
                    features_chrom[line[0][3:]].append([int(line[1]), int(line[2]), line[3]])
                except:
                    features_chrom[line[0][3:]].append([int(line[1]), int(line[2]) ])
                         
    return(features_chrom)


def TSS_matrix(annotation, tbx, barcodes):
    """
    Returns a scipy.sparse.lil_matrix with dimensions: nb_barcodes x nb_TSS 
    """
    nb_barcodes = len(barcodes)
    dict_barcodes = {barcodes[i]: i for i in range(0,len(barcodes))}

    window_list = []

    for chrom in sorted(annotation.keys()):
        window_list += [["".join(['chr', chrom]), int(n[0]), int(n[1])] for n in annotation[chrom]]

        mtx = lil_matrix((nb_barcodes, len(window_list)), dtype=np.float32) ## shape is NB_BARCODES x NB_TSS
        for i, tmp_feat in enumerate(window_list):
            for row in tbx.fetch(tmp_feat[0], tmp_feat[1], tmp_feat[2], parser=pysam.asTuple()):
                bcd = str(row).split('\t')[-2]
                if bcd not in barcodes:
                    continue
                else:
                    mtx[dict_barcodes[bcd], i] += 1

    return mtx



SLOP_TSS = "slop_tss.bed"
LEFT_FLANK = "lflank_tss.bed"
RIGHT_FLANK = "rflank_tss.bed"
META = "/proj/tmp/mukund/atac/AGG_ATAC_210420__NO_NORM/outs/singlecell.csv"
FRAGMENTS = "/proj/tmp/mukund/atac/AGG_ATAC_210420__NO_NORM/outs/fragments.tsv.gz"

tss_feat = load_features(SLOP_TSS)
lflank_feat = load_features(LEFT_FLANK)
rflank_feat = load_features(RIGHT_FLANK)

tbx = pysam.TabixFile(FRAGMENTS)

df = pd.read_csv(META)
df_filtered = df[(df.is__cell_barcode == 1)]
barcodes = list(set(df_filtered.barcode.tolist()))

with open("data/Tss_files/barcodes.pickle", "wb") as handle:
    pickle.dump(barcodes, handle, protocol=pickle.HIGHEST_PROTOCOL)

tss_mtx = TSS_matrix(tss_feat, tbx, barcodes).sum(axis=1)
with open("data/Tss_files/Agg_tss_mtx.pickle", "wb") as handle:
    pickle.dump(tss_mtx, handle, protocol=pickle.HIGHEST_PROTOCOL)

lflank_mtx = TSS_matrix(lflank_feat, tbx, barcodes).sum(axis=1)
with open("data/Tss_files/Agg_lflank_mtx.pickle", "wb") as handle:
    pickle.dump(lflank_mtx, handle, protocol=pickle.HIGHEST_PROTOCOL)

rflank_mtx = TSS_matrix(rflank_feat, tbx, barcodes).sum(axis=1)
with open("data/Tss_files/Agg_rflank_mtx.pickle", "wb") as handle:
    pickle.dump(rflank_mtx, handle, protocol=pickle.HIGHEST_PROTOCOL)

avg_noise = (lflank_mtx + rflank_mtx)/2
tss_score = (tss_mtx/avg_noise)

TSS_scores = [x[0] for x in tss_score.tolist()]

with open("data/Tss_files/Agg_TSS_scores.pickle", "wb") as handle:
    pickle.dump(TSS_scores, handle, protocol=pickle.HIGHEST_PROTOCOL)

