#!/bin/python

import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

## Check the correlation (relative distance) between peaks form different modalities

base = "/data/proj/GCB_MK/scCT/nanoCT_EAE/fragments/P27505_OLG"
os.system(f"chdir {base}")

## check relative distance vs frequency distribution
os.system("bedtools reldist -a H3K27ac.peaks -b ATAC.peaks > H3K27ac_v_ATAC_peaks_relDist.tsv ")
os.system("bedtools reldist -a H3K27ac.peaks -b H3K27me3.peaks > H3K27ac_v_H3K27me3_peaks_relDist.tsv")
os.system("bedtools reldist -a ATAC.peaks -b H3K27me3.peaks > ATAC_v_H3K27me3_peaks_relDist.tsv")
os.system("bedtools reldist -a peaks/human_OLIG2nuclei_H3K27ac_epilepsy_pooled_hg38.ucsc.peak -b H3K27me3.peaks > Nott_v_H3K27me3_peaks_relDist.tsv")
os.system("bedtools reldist -a peaks/human_OLIG2nuclei_H3K27ac_epilepsy_pooled_hg38.ucsc.peak -b H3K27ac.peaks > Nott_v_H3K27ac_peaks_relDist.tsv")


## plot summary
df = pd.read_csv(f"{base}/H3K27ac_v_ATAC_peaks_relDist.tsv", sep="\t")
df["mat"] = "H3K27ac_v_ATAC"
ax = sns.lineplot(data=df, x="reldist",y="fraction", alpha=0.8, hue="mat", palette=["red"], style="mat", lw=2);

df = pd.read_csv(f"{base}/H3K27ac_v_H3K27me3_peaks_relDist.tsv", sep="\t")
df["mat"] = "H3K27me3_v_H3K27ac"
ax = sns.lineplot(data=df, x="reldist",y="fraction", alpha=0.8, hue="mat", palette=["blue"], lw=2)

df = pd.read_csv(f"{base}/ATAC_v_H3K27me3_peaks_relDist.tsv", sep="\t")
df["mat"] = "ATAC_v_H3K27me3"
sns.lineplot(data=df, x="reldist",y="fraction", alpha=0.8, hue="mat", palette=["green"], lw=2)

df = pd.read_csv(f"{base}/Nott_v_H3K27ac_peaks_relDist.tsv", sep="\t")
df["mat"] = "Nott_et_al_v_H3K27ac"
sns.lineplot(data=df, x="reldist",y="fraction", alpha=0.8, hue="mat", palette=["orange"], lw=2)

df = pd.read_csv(f"{base}/Nott_v_H3K27me3_peaks_relDist.tsv", sep="\t")
df["mat"] = "Nott_et_al_v_H3K27me3"
sns.lineplot(data=df, x="reldist",y="fraction", alpha=0.8, hue="mat", palette=["pink"], lw=2)


for i in range(len(ax.lines)):
    ax.lines[i].set_linestyle("-.")
    
ax.legend(fontsize=10, frameon=False)
plt.suptitle('Relative distance between peak sets', fontsize=15)
plt.xlabel('Relative Distance', fontsize=12)
plt.ylabel('Frequency', fontsize=12)
ax.tick_params(axis='both', which='major', labelsize=10)
sns.despine(right=True, top=True)

plt.savefig(f"{base}/relDist_between_peak_sets.pdf")