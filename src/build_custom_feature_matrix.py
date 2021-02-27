#!/bin/python

import os
import sys
import scanpy as sc
import anndata as ad
import numpy as np
import pandas as pd
import episcanpy.api as epi
import argparse

## episcanpy environment should be activated

parser = argparse.ArgumentParser(description="Build custom feature matrix")
parser.add_argument("-f", "--fragments", help="path to fragments.tsv.gz file")
parser.add_argument("-p", "--peaks", help="path to the peak/features BED file")
parser.add_argument("-c", "--csv", help="path to the FORMATTED singlecell.csv file.")
parser.add_argument("-o", "--out", help="name of the output .h5ad file containing the matrix")
args = parser.parse_args()

## check conda environment
if sys.executable.split('/')[-3] == "episcanpy":
    pass
else:
    print("activate conda environment\n")
    sys.exit(1)
    
#check arguments
if not os.path.exists(args.fragments):
    print("fragments file not in the specified path\n")
    sys.exit(1)
if not os.path.exists(args.peaks):
    print("peaks file not in the specified path\n")
    sys.exit(1)
if not os.path.exists(args.csv):
    print("singlecell.csv file not in the specified path\n")
    sys.exit(1)
    
# load features
annot = epi.ct.load_features(args.peaks)

# build the count matrix. Takes a lot of time. 
counts = epi.ct.bld_mtx_fly(tsv_file = args.fragments,
                            csv_file = args.csv,
                            annotation = annot,
                            save = args.out)