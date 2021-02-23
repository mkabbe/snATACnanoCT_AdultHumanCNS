#!/bin/python

import os 
import pandas as pd


CR_META = pd.read_csv("data/P20057_1004/singlecell.csv", sep=",")
CR_META = CR_META[CR_META.is__cell_barcode == 1]
CR_BCD = CR_META.barcode.to_list()

MK_META = pd.read_csv("data/P20057_1004/MK_pre_singlecell.csv")
MK_BCD = MK_META.barcode.to_list()

print("Cells from cellranger: {} cells.".format(len(CR_BCD)))
print("Cells from MK_pre-processing: {} cells.".format(len(MK_BCD)))
print("Common cells: {} cells.".format(len([i for i in CR_BCD if i in MK_BCD])))
print("Unique to cellranger: {} cells.".format(len([i for i in CR_BCD if i not in MK_BCD])))
print("Unique to MK_pre-processing: {} cells.".format(len([i for i in MK_BCD if i not in CR_BCD])))


