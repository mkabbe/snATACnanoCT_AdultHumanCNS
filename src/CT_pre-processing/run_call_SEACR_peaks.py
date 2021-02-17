#!/bin/python

import os

samples = ["P20057_1004", "P20057_1005", "P20057_1006"]

for sample in samples:
        if not os.path.exists(sample):
                os.mkdir(sample)
        RUN_SEACR = "sbatch scripts/call_SEACR_peaks.sh {}".format(sample)
        os.system(RUN_SEACR)