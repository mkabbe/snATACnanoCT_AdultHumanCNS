import os
import anndata as ad
import scanpy as sc
import pandas as pd
import numpy as np
from scipy.stats import ks_2samp
import matplotlib.pyplot as plt
import seaborn as sns
#%config InlineBackend.figure_format = "retina"


## Plot signal distribution and run Kolmogorov-Smirnov test for significance

def smooth_array(x,y,s=5):    
    x_sm = np.array(x)
    y_sm = np.array(y,np.float64) #.replace(np.nan, 0))
    x_smooth = np.linspace(x_sm.min(), x_sm.max(), 200)
    spl = interpolate.UnivariateSpline(x, y.replace(np.nan, 0))
    sigma = s
    x_g1d = ndimage.gaussian_filter1d(x_sm, sigma)
    y_g1d = ndimage.gaussian_filter1d(y_sm, sigma)
    return (x_g1d, y_g1d)
    
    
hox = "HOXA"
    
df = pd.read_csv(f"/date/gcb/GCB_MK/profilePlots/f{hox}_50kb_OLG.profile_data.tsv", sep="\t").T
df.index = [x for x in range(-1,df.shape[0]-1)]
df = df.drop(index=[0]).fillna(0)
df = df.drop(index=[x for x in range(101,df.shape[0])])
#df.fillna(0)

df_atac = df[[1]].drop(index=[-1])
df_ac = df[[3]].drop(index=[-1])
df_me = df[[5]].drop(index=[-1])
df_ac.columns = [x for x in range(1,2)]
df_me.columns = [x for x in range(1,2)]


g = 10  #bin_size


fig, ax = plt.subplots(nrows=nrows,ncols=ncols, figsize=[6,3]) 

sig = 3 ## Smoothing to apply (for visualization ONLY - all stats are calculated on Raw data)
alpha = 0.05 # pval threshold
sig_windows = []
sig_values = []

col_dict = {0:"blue",1:"green",2:"red"}
for i,mod in enumerate([df_atac, df_ac, df_me]):

    xraw, yraw = mod.index, mod[1]

    for x in range(0,len(yraw),g):
        y1 = yraw[x:x+g]
        y2 = yraw[x+g:x+g+g]
        # Run the Kolmogorov-Smirnov test
        try:
            ks_stat, p_value = ks_2samp(y1, y2)
        except ValueError:
            pass

        if p_value < alpha:
            sig_windows.append(x+g) # The window at the border of the differing signal
            sig_values.append(p_value)

    xg1d, yg1d = smooth_array(xraw, yraw, sig)


    sns.lineplot(x=xg1d, y=yg1d, linewidth=1.5, color=col_dict[i], alpha=1, ax=ax)
    

for x in range(0,len(yg1d),2):
    plt.axvline(x, color="gray", lw=0.1,alpha=0.5)
    
for x in sig_windows:
    plt.axvline(x, color="black", lw=17, alpha=0.09)


ax.set_xlim([1,99])
ax.set_xticks([])

sns.despine()
plt.savefig(f"/date/gcb/GCB_MK/profilePlots/plots/{hox}_coverage_KSTest.png", dpi=400, bbox_inches="tight")

