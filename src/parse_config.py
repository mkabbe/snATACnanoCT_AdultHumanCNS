#!/bin/python

import os
import json 


with open("cr.json") as f:
    data = json.load(f)


with open("scATAC_AGGR.csv", "w") as in_:
    in_.write("library_id,fragments,cells,peaks\n")
    for sample in data["samples"]:
        file_prefix = "{0}/{1}/outs".format(data["cr_count_dir"],sample)
        in_.write("{0},{1}/fragments.tsv.gz,{1}/singlecell.csv,{1}/peaks.bed\n".format(sample, file_prefix))

