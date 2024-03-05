import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.colors import Normalize
from scipy.interpolate import interpn


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

    ax.scatter( x, y, c=z, s=6, **kwargs )

    norm = Normalize(vmin = np.min(z), vmax = np.max(z))
 #   cbar = fig.colorbar(cm.ScalarMappable(norm = norm), ax=ax)
 #   cbar.ax.set_ylabel('Density')

    return ax



OUT_PATH = "/data/proj/GCB_MK/CZI/episcanpy_analysis/AGG_ATAC_210218/data/SAMPLE_DATA/"

samples = pd.read_csv(OUT_PATH + "QC_stats_NEW.csv")
samples = list(samples["NGI_ID"])

fig, axes = plt.subplots(6, 8, figsize=(16, 12))

axes = axes.flatten()

for i, ax in enumerate(axes):
    if i < len(samples):
        df = pd.read_csv(OUT_PATH+ samples[i] + "_TSS_obs.csv")
        df["TSSe_thresh"] = [min(30,x) for x in df["TSS_score"]]
        ax = density_scatter(x = df.log_nb_features, y=df.TSSe_thresh, ax = ax)
        ax.set_xlabel("log_nb_features")
        ax.set_ylabel("TSS enrichment")
        ax.set_title(samples[i])

# Adjust spacing between subplots
fig.tight_layout()
#Save
plt.savefig(OUT_PATH+ "all_TSSvFrag_QC_plot.png", dpi=300, bbox_inches="tight")