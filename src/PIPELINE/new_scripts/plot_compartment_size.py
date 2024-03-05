#!/bin/python


import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats



df = pd.read_csv("/date/gcb/GCB_MK/microC/bam/hOPC_100kb_pca1.bedgraph", sep="\t", header=None)
df = df.rename(columns={0:"chrom", 1:"start", 2:"end", 3:"score"})
df["region"] = df.chrom.astype(str) + "_" + df.start.astype("str") + "_" + df.end.astype(str)
df["AB"] = ["A" if x>0 else "B" for x in df.score]
region_dict = dict(zip(df.region, df.AB))
chroms = list(set(df.chrom))
regions = list(region_dict.keys())


ab_border = []
for chrom in chroms:
    for i in range(len(regions)):
    #for i,region in enumerate(regions):
        if regions[i].split("_")[0] != chrom:
            continue
        else:
            if i==0:
                ab_border.append(region)
                continue
            try:
                if region_dict[regions[i]] != region_dict[regions[i+1]]:
                    ab_border.append((regions[i],region_dict[regions[i]]))
            except IndexError:
                ab_border.append((regions[i],region_dict[regions[i]])) #last bin in chromosome
                
df["is_border"] = False
df.loc[df.region.isin(set([x[0] for x in ab_border])),"is_border"] = True


plt.rcParams["figure.figsize"] = [4,6]
fig = plt.figure()
ax = fig.add_subplot(2, 1, 1)

comp_size = []
for chrom in set(df.chrom):

    df_border=df.loc[(df.is_border==True)&(df.chrom==chrom)]
#    comp_size = []
    for idx in range(len(df_border["end"])):
        if idx==0:
            comp_size.append(0)
        else:
            comp_size.append(df_border["end"].tolist()[idx]-df_border["end"].tolist()[idx-1])

count, bins_count = np.histogram(comp_size, bins=500)
pdf = count / sum(count)
cdf = np.cumsum(pdf)
# plotting CDF
ax.plot(bins_count[1:], cdf, linewidth=2, alpha=0.8,color="black")
#ax.set_xscale('log')
#sns.despine()
#sns.displot(comp_size,kind="kde",fill=True, alpha=0.1)

#ax.set_xlim(0, 70)
ax.set_xticks([0,6.2e6,2e7,4e7,6e7,8e7])
ax.set_xticklabels(['0','*','20','40','60','80'])

ax.set_xlabel("Compartment size (in Mb)");
ax.set_ylabel("Cumulative fraction");

ax.axhline(y=0.95, linestyle="--", linewidth=1, color="silver")
ax.axvline(x=6.2e6, linestyle="--", linewidth=1, color="silver")