import os
import csv
import scipy
import scanpy as sc
import anndata as ad
import numpy as np
import pandas as pd
import episcanpy.api as epi
from sklearn.feature_extraction.text import TfidfTransformer
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib import style
import warnings
import magic
style.use("default")
warnings.filterwarnings('ignore')
%config InlineBackend.figure_format = "retina"
os.chdir("/data/proj/GCB_MK/CZI/episcanpy_analysis/AGG_ATAC_210218")



ol_adata = ad.read("data/feature_matrices/210617_OLG_LSI-2_results_matrix.h5ad")



## split ages into different bins, not needed for this specific analysis
#ol_adata.obs['age_bins'] = pd.cut(x=ol_adata.obs['Age'],
#                right=False, bins=[25, 35, 45, 55, 65, 75],
#                labels=["25-34","35-44","45-54","55-64","65-74"])



## Compare OPCs from CSC and BA4
ol_subdata = ol_adata[(ol_adata.obs["annot_rough"]=="OPC_1") & (ol_adata.obs["Tissue"]!="CB")]


epi.tl.rank_features(ol_subdata, groupby="Tissue")
subresult = ol_subdata.uns["rank_features_groups"]
subgroups = subresult["names"].dtype.names

opc_csc_markers = subresult["names"]["CSC"][:10].tolist()
opc_ba_markers = subresult["names"]["BA4"][:10].tolist()

####     plots Trackplot in jupyter
#ax = sc.pl.tracksplot(ol_subdata, opc_csc_markers, groupby='Tissue')
#ax = sc.pl.tracksplot(ol_subdata, opc_ba_markers, groupby='Tissue')

####     plots Heatmap in jupyter
#sc.pl.rank_genes_groups_heatmap(ol_subdata, key="rank_features_groups", vmax="4");

opc_csc_var_peaks = pd.DataFrame(
        {key: subresult[key]["CSC"] 
         for key in ["names", "logfoldchanges", "pvals", "pvals_adj"]})

opc_ba_var_peaks = pd.DataFrame(
        {key: subresult[key]["BA4"] 
         for key in ["names", "logfoldchanges", "pvals", "pvals_adj"]})


opc_csc_var_peaks_logfold = opc_csc_var_peaks.sort_values(by=["logfoldchanges"], ascending=False).head(20).names.tolist()
#opc_csc_var_peaks.sort_values(by=["logfoldchanges"], ascending=False).head(20)

opc_ba_var_peaks_logfold = opc_ba_var_peaks.sort_values(by=["logfoldchanges"], ascending=False).head(20).names.tolist()
#opc_csc_var_peaks.sort_values(by=["logfoldchanges"], ascending=False).head(20)


####     plots Trackplot in jupyter
#ax = sc.pl.tracksplot(ol_subdata, opc_csc_var_peaks_logfold, groupby='Tissue')
#ax = sc.pl.tracksplot(ol_subdata, opc_ba_var_peaks_logfold, groupby='Tissue')

csc_var_peak_df = opc_csc_var_peaks.sort_values(by=["logfoldchanges"], ascending=False)ba_var_peak_df = opc_ba_var_peaks.sort_values(by=["logfoldchanges"], ascending=False)
ba_var_peak_df["logFC"] = [-x for x in ba_var_peak_df["logfoldchanges"]]
del ba_var_peak_df["logfoldchanges"]
ba_var_peak_df = ba_var_peak_df.rename(columns={"logFC":"logfoldchanges"})


df = pd.concat([csc_var_peak_df, ba_var_peak_df])
df["-log10(pval_adj)"] = [-1*np.log10(x) for x in df.pvals_adj]

df["significant"] = [False for _ in range(len(df))]
df.loc[(df["-log10(pval_adj)"] >= 3) & ((df["logfoldchanges"] >= 1.5) | (df["logfoldchanges"] <= -1.5)), "significant"] = True
df = df.rename(columns={"logfoldchanges":"logFC"})


####    Plots Volcanoplot in jupyter
#plt.figure(figsize=(8,5))
#ax= sns.scatterplot(data=df, x="logFC", y="-log10(pval_adj)", s=9, hue="significant", palette=["silver", "xkcd:violet"]);
#plt.axhline(y=3, alpha=0.1, ls="--", color='black', linewidth=0.9); ## pval_adj cutoff = 0.001
#plt.axvline(x=1.5, alpha=0.1, ls="--", color='black', linewidth=0.9);
#plt.axvline(x=-1.5, alpha=0.1, ls="--", color='black', linewidth=0.9);
#ax.get_legend().remove()


sig_df = df.loc[df["significant"]==True]
sig_df = sig_df.sort_values(by=["logFC"], ascending=False)


df_1 = sig_df.names.str.split("_", expand=True)
df_1["name"] = ["peak"+str(x) for x in range(len(df_1))]
df_1["strand"] = ["+" for x in range(len(df_1))]
df_1["logFC"] = sig_df["logFC"]
df_1["pvals"] = sig_df["pvals"]
df_1["pvals_adj"] = sig_df["pvals"]

df_1.to_csv("data/variable_peaks/sig_diff_peaks_OPC_CSCvsBA4.bed", header=False, index=False, sep="\t")

### BED FILE THAT WAS WRITTEN IS THE INPUT FOR HOMER ANNOTATION (annotate_peaks_homer.sh)
