#!/bin/python

import anndata as ad
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
%config InlineBackend.figure_format = "retina"

adata = ad.read("data/epi/top60k_binned_2kb_count_matrix_RESULTS.h5ad")

def plot_reads_by_tissue(adata, out_file):
    adata = ad.read("data/epi/top60k_binned_2kb_count_matrix_RESULTS.h5ad")
    tissue_types = ["CB", "CSC", "BA4"]  # this is invariate so can be hard-coded
    feat_dct = {}
    
    for tissue in tissue_types:
        feat_dct[tissue] = adata.obs.loc[adata.obs.Tissue == tissue, "log_nb_features"]

    CSC = feat_dct["CSC"]
    CB = feat_dct["CB"]
    BA4 = feat_dct["BA4"]

    fig, ax = plt.subplots(1, 1,figsize=(6,6))
    sns.boxplot(x='Tissue',y='log_nb_features',
                data=pd.DataFrame({'log_nb_features':pd.concat([CSC,BA4,CB],axis=0),
                                   'Tissue':np.repeat(["CSC","BA4","CB"],[len(CSC),len(BA4),len(CB)])}),
                ax = ax,palette=colors,width=0.5)
    plt.savefig(out_file)




