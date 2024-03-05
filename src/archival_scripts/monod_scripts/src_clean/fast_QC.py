import os
import sys
import scanpy as sc
import anndata as ad
import numpy as np
import pandas as pd
import episcanpy.api as epi
import pysam
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.colors import Normalize
from scipy.interpolate import interpn
from scipy.sparse import lil_matrix

## episcanpy environment should be activated
"""
Builds a 2kb count matrix, gets TSS_score and generates a TSS_score_vs_log_nb_features density scatter plot

TODO : 
-> automate filtering out of cells with <1000 features and TSS_score < 3. 
-> Update the singlecell.csv[is__cell_barcode] attribute.
-> Note samples where majority of cells have poor QC metrics. These should be removed before running cellranger-atac aggr. 

"""
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
    set_barcode = set(barcodes)

    for chrom in sorted(annotation.keys()):
        window_list += [["".join(['chr', chrom]), int(n[0]), int(n[1])] for n in annotation[chrom]]

        mtx = lil_matrix((nb_barcodes, len(window_list)), dtype=np.float32) ## shape is NB_BARCODES x NB_TSS
        for i, tmp_feat in enumerate(window_list):
            for row in tbx.fetch(tmp_feat[0], tmp_feat[1], tmp_feat[2], parser=pysam.asTuple()):
                bcd = str(row).split('\t')[-2]
                if bcd not in set_barcode:
                    continue
                else:
                    mtx[dict_barcodes[bcd], i] += 1
    return mtx

def density_scatter( x , y, ax = None, sort = True, bins = 20, **kwargs )   :
    """
    Scatter plot of log_nb_fragments vs TSS_score (density colored)
    """
    if ax is None :
        fig , ax = plt.subplots()
    data , x_e, y_e = np.histogram2d( x, y, bins = bins, density = True )
    z = interpn( ( 0.5*(x_e[1:] + x_e[:-1]) , 0.5*(y_e[1:]+y_e[:-1]) ) , data , np.vstack([x,y]).T , method = "splinef2d", bounds_error = False)

    #To be sure to plot all data
    z[np.where(np.isnan(z))] = 0.0

    # Sort the points by density, so that the densest points are plotted last
    if sort :
        idx = z.argsort()
        x, y, z = x[idx], y[idx], z[idx]

    ax.scatter( x, y, c=z, **kwargs )

    norm = Normalize(vmin = np.min(z), vmax = np.max(z))
    cbar = fig.colorbar(cm.ScalarMappable(norm = norm), ax=ax)
    cbar.ax.set_ylabel('Density')

    return ax


DATA_PATH = "/data/proj/GCB_MK/10XATAC/data/scatac/cellranger_COUNT_output/"
#SAMPLE_PATH = DATA_PATH + sys.argv[1] + "/outs/"
#SAMPLE_PATH = sys.argv[1]
OUT_PATH = "/data/proj/GCB_MK/CZI/episcanpy_analysis/AGG_ATAC_210218/data/SAMPLE_DATA/"

FEATURES = "/data/proj/GCB_MK/CZI/episcanpy_analysis/AGG_ATAC_210218/ref/2kb_binned_hg38.bed"
#FRAGMENTS = SAMPLE_PATH + "fragments.tsv.gz"
FRAGMENTS = sys.argv[1]
CSV = FRAGMENTS.split("fragments")[0] + "singlecell_EPI.csv"
sample_id = FRAGMENTS.split("/")[-3] 
SLOP_TSS = "/data/proj/GCB_MK/CZI/episcanpy_analysis/AGG_ATAC_210218/ref/slop_tss.bed"
LEFT_FLANK = "/data/proj/GCB_MK/CZI/episcanpy_analysis/AGG_ATAC_210218/ref/lflank_tss.bed"
RIGHT_FLANK = "/data/proj/GCB_MK/CZI/episcanpy_analysis/AGG_ATAC_210218/ref/rflank_tss.bed"

# load features
annot = epi.ct.load_features(FEATURES)

# build the count matrix. Takes a lot of time.
adata = epi.ct.bld_mtx_fly(tsv_file = FRAGMENTS,
                            csv_file = CSV,
                            annotation = annot)

# load info to adata.obs 
epi.pp.filter_cells(adata, min_features=1)
adata.obs["log_nb_features"] = [np.log10(x) for x in adata.obs["nb_features"]]
adata.obs["barcode"] = adata.obs_names


tss_feat = load_features(SLOP_TSS) ## 50bp window
lflank_feat = load_features(LEFT_FLANK) ## 100bp window 2kb upstream
rflank_feat = load_features(RIGHT_FLANK) ## 100bp window 2kb downstream

tbx = pysam.TabixFile(FRAGMENTS)

barcodes = adata.obs["barcode"].tolist()

# build the TSS matrix
tss_mtx = TSS_matrix(tss_feat, tbx, barcodes).sum(axis=1)
lflank_mtx = TSS_matrix(lflank_feat, tbx, barcodes).sum(axis=1)
rflank_mtx = TSS_matrix(rflank_feat, tbx, barcodes).sum(axis=1)

avg_noise = (lflank_mtx + rflank_mtx)/2

for i in avg_noise:
    i[0] = max(i[0],0.2) ## lower threshold of 0.2 for low signal samples

score = tss_mtx/avg_noise

## upper threshold of 30
adata.obs["TSS_score"] = [min(round(i[0],3),30) for i in score.tolist()]

## write to csv file 
adata.obs.to_csv(OUT_PATH+ sample_id + "_TSS_obs.csv", header=True, index_label=False)

## plot and save figure
fig = plt.figure()
ax1 = density_scatter(x = adata.obs.log_nb_features, y = adata.obs.TSS_score)
ax1.set_xlabel("log_nb_features")
ax1.set_ylabel("TSS enrichment")
plt.savefig(OUT_PATH+ sample_id + "_TSSvFrag_QC_plot.pdf")
