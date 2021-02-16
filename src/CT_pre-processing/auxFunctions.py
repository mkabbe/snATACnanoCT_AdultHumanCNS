#!/usr/bin/python

import os 
import gzip
from shutil import move 


# remove whitespace and reformat for downstream applications

def stripReformat(in_file, TMP="tmp.txt"):
    with open(in_file, "r") as _in:
        with open(TMP, "a+") as _TMP:
            for line in _in:
                _TMP.write("\t".join(line.strip().split(" "))+"\n")
    move(TMP, in_file)
    return in_file



# filter fragments file to only contain fragments from cells passing CT filters

def filterFragments(barcodelist, fragments, out_gz):
    with gzip.open(fragments, "rt") as handle, open(out_gz, "a+") as _out:
        for line in handle:
            if any(barcode in line for barcode in barcodelist):
                _out.write(line+"\n")
