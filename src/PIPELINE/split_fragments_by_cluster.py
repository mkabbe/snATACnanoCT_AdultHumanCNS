## adapted from MB's scCT pipeline
import os
import anndata as ad
import pandas as pd
import pysam
import argparse
import gzip
import csv


parser = argparse.ArgumentParser(description="Extract Fragments for cells in each cluster")
parser.add_argument("-a", "--anndata", help="anndata file")
parser.add_argument("-f", "--fragments", action="store")
parser.add_argument("-b", "--barcodes", action="store")
parser.add_argument("-o", "--output", action="store")
args = parser.parse_args()


def parse_barcodes(args):
    bcd = {}
    with open(args.barcodes) as bc_csv:
        barcodes = csv.reader(bc_csv)
        for row in barcodes:
            bcd[row[0]] = row[1]
    return bcd


def filter_fragments_by_cluster(args, bcd):
    cluster_frag = {}
    with gzip.open(args.fragments) as f:
        tbx = pysam.tabix_iterator(f,pysam.asBed())
        for line in tbx:
            if line.name in bcd:
                try:
                    cluster_frag[bcd[line.name]].append(str(line))
                except KeyError:
                    cluster_frag[bcd[line.name]] = [str(line)]
    return cluster_frag


def main(args):
    
    adata = ad.read(args.anndata)
    if "barcode" not in adata.obs.columns:
        adata.obs["barcode"] = adata.obs.index.copy()

    try:
        bc_cluster_df = adata.obs[["barcode","leiden"]]
    except KeyError:
        bc_cluster_df = adata.obs[["barcode","louvain"]]
    bc_cluster_df.to_csv(args.barcodes, header=False, index=False)
    
    barcodes  = parse_barcodes(args)
    reads_dic = filter_fragments_by_cluster(args.barcodes)
    for key in reads_dict:
        with open(args.output + key + ".bed",'w') as f:
            f.write("\n".join(reads_dic[key]))

main(args)