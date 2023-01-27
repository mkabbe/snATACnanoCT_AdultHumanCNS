#!/bin/python

import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import os
import cooltools
import cooler
from matplotlib.ticker import EngFormatter
from matplotlib.colors import LogNorm
import cooltools.lib.plotting
import itertools

#mcooler_file = "P22959_hOPC.contact-map.mcool"

cooler_file = "P27454_REdeep_hOPC.mapped.pairs.contact-map.100000.cool"
hicooler_file = "P27454_REdeep_hOPC.mapped.pairs.contact-map.25000.cool"
supercooler_file = "P27454_hOPC.mapped.pairs.contact-map.5000.cool"
res = 1000000
hi_res = 25000
super_res = 5000

#cooler.fileops.list_coolers(mcooler_file)
#hopc_clr = cooler.Cooler(f"{mcooler_file}::resolutions/{res}")
hopc_clr = cooler.Cooler(f"{cooler_file}")

chromstarts = []
for i in hopc_clr.chromnames:
    print(f'{i} : {hopc_clr.extent(i)}')
    chromstarts.append(hopc_clr.extent(i)[0])

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

f, axs = plt.subplots(
    figsize=(14,4),
    ncols=3)

#ax = axs[0]
#im = ax.matshow(hopc_clr.matrix(balance=False)[:], vmax=30, cmap="fall");
#plt.colorbar(im, ax=ax ,fraction=0.046, pad=0.04, label='raw counts');
#ax.set_xticks(chromstarts)
#ax.set_yticks(chromstarts)
#ax.xaxis.tick_bottom()
#ax.set_title('All chromosomes')



### PLOT whole Chromosome
ax = axs[0]



im = ax.matshow(hopc_clr.matrix(balance=False).fetch("chr7"), vmax=200,cmap="bwr",
                extent=(0,hopc_clr.chromsizes['chr7'],
                hopc_clr.chromsizes['chr7'], 0));

plt.colorbar(im, ax=ax ,fraction=0.046, pad=0.04, label='raw counts');

ax.set_title('chr7', y=1.08)
ax.set_ylabel('position, Mb')
format_ticks(ax)

### PLOT Part of Chromosome
ax = axs[1]
start, end = 25_000_000, 30_000_000
region = ('chr7', start, end)

#supres_hopc_clr = cooler.Cooler(f"{mcooler_file}::resolutions/{super_res}")
hires_hopc_clr = cooler.Cooler(f"{hicooler_file}"})
im = ax.matshow(
    hires_hopc_clr.matrix(balance=False).fetch(region), vmax=40, cmap="bwr",
    extent=(start, end, end, start)
);
ax.set_title(f'chr7:{start:,}-{end:,}', y=1.08)
plt.colorbar(im, ax=ax ,fraction=0.046, pad=0.04, label='raw counts');
format_ticks(ax)

### PLOT smaller part of Chromosome
ax = axs[2]
start, end = 26_500_000, 28_500_000
region = ('chr7', start, end)

#super_hopc_clr = cooler.Cooler(f"{mcooler_file}::resolutions/{super_res}")
super_hopc_clr = cooler.Cooler(f"{supercooler_file}")

norm = LogNorm(vmax=0.01, vmin=0.0001) ## set upper and lower limits of color scale
data = super_hopc_clr.matrix(balance=True).fetch(region)
im = pcolormesh_45deg(ax, data, start=region[1], resolution=super_res, norm=norm, cmap='bwr')

#im = ax.matshow(
#    super_hopc_clr.matrix(balance=False).fetch(region), cmap="bwr", vmax=30,
#    extent=(start, end, end, start)
#);

ax.set_title(f'chr7:{start:,}-{end:,}', y=1.08)
plt.colorbar(im, ax=ax ,fraction=0.046, pad=0.04, label='raw counts');
format_ticks(ax)



plt.tight_layout()
#plt.savefig("hOPC_HOXA_cluster_NEW.pdf")
#plt.savefig("hOPC_HOXA_cluster_TRIAL.png")
plt.savefig("bCell_HOXA_cluster_TRIAL.png")
