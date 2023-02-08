#!/bin/python
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import os
import sys
import cooltools
import cooler
from matplotlib.ticker import EngFormatter
from matplotlib.colors import LogNorm
from mpl_toolkits.axes_grid1 import make_axes_locatable
import cooltools.lib.plotting
import itertools
import argparse
from matplotlib.colors import LinearSegmentedColormap



coolwarm_extremes = LinearSegmentedColormap.from_list('coolwarm_extremes', (
        # Edit this gradient at https://eltos.github.io/gradient/#coolwarm_extremes=0:400083-3:7414D5-4.8:7038DB-10.8:1201FF-14.9:2D9ADA-24.6:1DCDEE-40:B3CAD8-57.3:B3CAD8-68.1:D4B7DC-76:DC71DA-87.5:FF0200-92.2:E47709-95:EB8D2A-100:FFFA00
        (0.000, (0.251, 0.000, 0.514)),
        (0.030, (0.455, 0.078, 0.835)),
        (0.048, (0.439, 0.220, 0.859)),
        (0.108, (0.071, 0.004, 1.000)),
        (0.149, (0.176, 0.604, 0.855)),
        (0.246, (0.114, 0.804, 0.933)),
        (0.400, (0.702, 0.792, 0.847)),
        (0.573, (0.702, 0.792, 0.847)),
        (0.681, (0.831, 0.718, 0.863)),
        (0.760, (0.863, 0.443, 0.855)),
        (0.875, (1.000, 0.008, 0.000)),
        (0.922, (0.894, 0.467, 0.035)),
        (0.950, (0.922, 0.553, 0.165)),
        (1.000, (1.000, 0.980, 0.000))))
# cmap = coolwarm_extremes

parser = argparse.ArgumentParser(add_help=True, description='Plot Micro-C contact matrix')
    
parser.add_argument("--cool", "-m", help=".cool matrix", required=True)
parser.add_argument("--res", "-r", help="resolution", required=True, type=int)
parser.add_argument("--chr", "-c", help="chromosome", required=True)
parser.add_argument("--start", "-s",  help="start position", required=True, type=float)
parser.add_argument("--end", "-e",  help="end position", required=True, type=float)
parser.add_argument("--output", "-o",  help="file name for plot", required=True)
parser.add_argument("--balanced", help="is matrix balanced (bool)", default=True)  
parser.add_argument("--vmax", help="max value for cmap", default=0.01, type=float)
parser.add_argument("--vmin", help="min value for cmap", default=0, type=float)
parser.add_argument("--cmap", help="colormap", default="RdYlBu")
parser.add_argument("--title", help="Plot title; default will be region")
parser.add_argument("--dpi", help="dpi", default=400, type=int)

args = parser.parse_args()


def pcolormesh_45deg(ax, matrix_c, start=0, resolution=1, *args, **kwargs):
    start_pos_vector = [start+resolution*i for i in range(len(matrix_c)+1)]
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
    

def format_ticks(ax, x=True, y=True, rotate=True):
    if y:
        ax.yaxis.set_major_formatter(bp_formatter)
    if x:
        ax.xaxis.set_major_formatter(bp_formatter)
        ax.xaxis.tick_bottom()
    if rotate:
        ax.tick_params(axis='x',rotation=45)


bp_formatter = EngFormatter('b')


super_hopc_clr = cooler.Cooler(args.cool)

chromstarts = []
for i in super_hopc_clr.chromnames:
    #print(f'{i} : {hopc_clr.extent(i)}')
    chromstarts.append(super_hopc_clr.extent(i)[0])

if args.dpi is not None:
    plt.rcParams['figure.dpi'] = int(args.dpi)
    
f, ax = plt.subplots(figsize=(4,4))

chrom = args.chr
start = args.start
end = args.end
region = (chrom, start, end)

if args.title is not None:
    title = args.title
else:
    title = f"{chrom}:{int(start)}-{int(end)}"

norm = LogNorm()    # vmax=0.1, vmin=0.01) ## set upper and lower limits of color scale
data = super_hopc_clr.matrix(balance=args.balanced).fetch(region)
im = pcolormesh_45deg(ax, data, start=region[1], resolution=args.res, vmax=args.vmax, vmin=args.vmin, cmap=args.cmap)
#im = pcolormesh_45deg(ax, data, start=region[1], resolution=args.res, cmap=args.cmap) ## for obs/exp
ax.set_aspect(0.5)
ax.set_ylim(0,end-start)

format_ticks(ax)
ax.set_title(title, y=1.08, fontsize=8)

sns.despine(bottom=True, left= True, trim=True, offset = 1)
plt.xticks(color='black', rotation=45, fontsize='3', horizontalalignment='right')
plt.tick_params(width=0.5, labelleft = False, left=False)


divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.05, aspect = 1000) ## contact-map colorbar position+size
cbar = plt.colorbar(im, cax=cax,)
cbar.set_label('ICE norm. counts',size=5)
cbar.ax.tick_params(labelsize=4)

plt.tight_layout()
plt.savefig(args.output)