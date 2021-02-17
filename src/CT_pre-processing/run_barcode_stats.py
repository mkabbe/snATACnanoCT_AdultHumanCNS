#!/bin/python

import os 

## for submitting barcode_statistics scripts 

samples = ["P20057_1004", "P20057_1005", "P20057_1006"]

if not os.path.exists("barcode_statistics"):
    os.mkdir("barcode_statistics")

for sample in samples:
    BCD_STATS_ALL = "sbatch scripts/barcode_stats_ALL.sh {}".format(sample)
    BCD_STATS_PEAK = "sbatch scripts/barcode_stats_PEAK.sh {}".format(sample)
    os.system(BCD_STATS_ALL)
    os.system(BCD_STATS_PEAK)
