#!/bin/python

import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import os
import cooltools
import cooler
from matplotlib.ticker import EngFormatter
import matplotlib.patches as patches
os.chdir("/data/proj/GCB_MK/micro-C")
bp_formatter = EngFormatter('b')
%config InlineBackend.figure_format = "retina"
import cooltools.lib.plotting


#mcooler_file = "P22959_hOPC.contact-map.mcool"
mcooler_file = "/date/gcb/GCB_MK/microC/bam/obs_exp_hOPC_100kb.cool"
#mcooler_file = "/date/gcb/GCB_MK/microC/bam/cools/hOPC_100kb.cool"

hires = 100000


hires_clr = cooler.Cooler(mcooler_file)

def format_ticks(ax, x=True, y=True, rotate=True):
    if y:
        ax.yaxis.set_major_formatter(bp_formatter)
    if x:
        ax.xaxis.set_major_formatter(bp_formatter)
        ax.xaxis.tick_bottom()
    if rotate:
        ax.tick_params(axis='x',rotation=45)

plt.rcParams['figure.dpi'] = 200

f, ax = plt.subplots(figsize=(4,4),ncols=1);
#ax = axs[55
im = ax.matshow(hires_clr.matrix(balance=False).fetch("chr7"), vmax=1, cmap="fall",
                extent=(0,hires_clr.chromsizes['chr7'],
                hires_clr.chromsizes['chr7'], 0));

plt.colorbar(im, ax=ax ,fraction=0.046, pad=0.04, label='ICE norm counts (obs/exp)');

ax.set_title('chr7', y=1.08)
ax.set_ylabel('position, Mb')
format_ticks(ax)

