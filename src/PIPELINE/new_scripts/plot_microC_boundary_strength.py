# import core packages
import os
import warnings
warnings.filterwarnings("ignore")
#from itertools import combination

# import semi-core packages
import matplotlib.pyplot as plt
from matplotlib import colors
%matplotlib inline
plt.style.use('seaborn-poster')
import numpy as np
import pandas as pd

# import open2c libraries
import bioframe

import cooler
import cooltools

from cooltools import insulation

resolution = 10000

clr_b = cooler.Cooler("cools/hOPC_10kb.cool")
windows = [3*resolution, 5*resolution, 10*resolution, 25*resolution]
insulation_table = insulation(clr_b, windows, verbose=True)

first_window_summary =insulation_table.columns[[ str(windows[-1]) in i for i in insulation_table.columns]]
#insulation_table[['chrom','start','end','region','is_bad_bin']+list(first_window_summary)].iloc[1000:1005]


# Functions to help with plotting
def pcolormesh_45deg(ax, matrix_c, start=0, resolution=1, *args, **kwargs):
    start_pos_vector = [start+resolution*i for i in range(len(matrix_c)+1)]
    import itertools
    n = matrix_c.shape[0]
    t = np.array([[1, 0.5], [-1, 0.5]])
    matrix_a = np.dot(np.array([(i[1], i[0])
                                for i in itertools.product(start_pos_vector[::-1],
                                                           start_pos_vector)]), t)
    x = matrix_a[:, 1].reshape(n + 1, n + 1)
    y = matrix_a[:, 0].reshape(n + 1, n + 1)
    im = ax.pcolormesh(x, y, np.flipud(matrix_c), *args, **kwargs)
    im.set_rasterized(True)
    return im

from matplotlib.ticker import EngFormatter
bp_formatter = EngFormatter('b')
def format_ticks(ax, x=True, y=True, rotate=True):
    if y:
        ax.yaxis.set_major_formatter(bp_formatter)
    if x:
        ax.xaxis.set_major_formatter(bp_formatter)
        ax.xaxis.tick_bottom()
    if rotate:
        ax.tick_params(axis='x',rotation=45)
        
        
from matplotlib.colors import LogNorm
from mpl_toolkits.axes_grid1 import make_axes_locatable
import bioframe
plt.rcParams['font.size'] = 3

#HOXA - 25mb-30mb chr7
#HOXD - 173.5mb-178.5mb chr2

start = 173_500_000
end = 178_500_000
region = ('chr2', start, end)
title = "chr2: 173.5 Mb - 178.5 Mb"

norm = LogNorm(vmax=0.01, vmin=0.0001) ## set upper and lower limits of color scale

data = clr_b.matrix(balance=True).fetch(region)
f, ax = plt.subplots(figsize=(18, 6))
im = pcolormesh_45deg(ax, data, start=region[1], resolution=resolution, norm=norm, cmap='bwr')
ax.set_aspect(0.5)
#ax.set_ylim(0, 10*windows[0])
#ax.set_ylim(0, end-start)#400_000)
ax.set_ylim(0,1_000_000)
format_ticks(ax, rotate=False)
ax.xaxis.set_visible(False)

divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="1%", pad=0.1, aspect=6) ## contact-map colorbar position+size
plt.colorbar(im, cax=cax)

insul_region = bioframe.select(insulation_table, region)

ins_ax = divider.append_axes("bottom", size="50%", pad=-0.4, sharex=ax) # position of 
ins_ax.set_prop_cycle(plt.cycler("color", plt.cm.plasma(np.linspace(0,1,5))))

## Show insulation plots for all windows
for res in windows[2:3]:
    ins_ax.plot(insul_region[['start', 'end']].mean(axis=1), insul_region[f'log2_insulation_score_{res}'], label=f'Window {res} bp', linewidth=2, color="grey")
ins_ax.legend(bbox_to_anchor=(0., -1), loc='lower left', ncol=4);

ax.set_title(title)

format_ticks(ins_ax, y=False, rotate=False)
ax.set_xlim(region[1], region[2])


boundaries = insul_region[~np.isnan(insul_region[f'boundary_strength_{windows[2]}'])]
weak_boundaries = boundaries[~boundaries[f'is_boundary_{windows[2]}']]
strong_boundaries = boundaries[boundaries[f'is_boundary_{windows[2]}']]
ins_ax.scatter(weak_boundaries[['start', 'end']].mean(axis=1),
            weak_boundaries[f'log2_insulation_score_{windows[2]}'], label='Weak boundaries', s=40,color="dodgerblue")
ins_ax.scatter(strong_boundaries[['start', 'end']].mean(axis=1),
            strong_boundaries[f'log2_insulation_score_{windows[2]}'], label='Strong boundaries', s=40, color="red")

ins_ax.legend(bbox_to_anchor=(0., -1), loc='lower left', ncol=4);

format_ticks(ins_ax, y=False, rotate=False)
ax.set_xlim(region[1], region[2])

f.savefig("/data/proj/GCB_MK/scCT/nanoCT_EAE/plots/microC_TAD_boundary_strength_chr2.png",dpi=400, bbox_inches="tight")

