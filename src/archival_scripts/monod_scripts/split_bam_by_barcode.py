#!/bin/python

## Creates small bam files for each barcode in the anndata object

import sys
import os
import pysam
#import anndata as ad   <- import inside the function. Not always needed globally
from contextlib import ExitStack


def barcodesToCsv(adata_path, outdir):
    """
    feeder function for splitBam()
    makes barcodes.csv files for each SAMPLE (NGI_ID code).
    barcodes.csv file will contain the aggregate Barcode and the cell barcode
    E.g P20056_1004_barcodes.csv --> AGCGGTT-4,AGCGGTT-1 ##line1
    """
    import anndata as ad
    adata = ad.read(adata_path)
    
    barcodes_dict = {}
    for _id_ in set(adata.obs["NGI_ID"]):
        barcodes_dict[_id_] = adata.obs.loc[adata.obs["NGI_ID"] == _id_, "barcode"].to_list()
    
    for _id in barcodes_dict.keys():
        with open("{0}/{1}_barcodes.csv".format(outdir,_id),"w") as f:
            for bcd in barcodes_dict[_id]:
                base,suffix = bcd.split("-")
                f.write(bcd + "," + bcd.split("-")[0] + "-1\n")
                


def splitBam(bamFile, barcode_annotations, out_prefix, barcode_prefix="NA"):
    """
    Refactored Marek's CUT&Tag code as a function.
    Makes bam files for each barcode
    bamFile: Possorted bam file from 10x with barcodes e.g. AAACGAAAGTTGCTTG-1
    barcode_annotations: cluster-barcode csv. e.g. mOL,H3K27ac_N1_AAACGAAAGTTGCTTG-1
    out_prefix: folder to write the output  e.g. /path/to/output
    barcode_prefix: Sample name to prefix to bam barcode e.g. mOL,H3K27ac_N1  (Use "NA" if no prefix is needed)
    """
    if barcode_prefix == "NA":
      barcode_prefix = ""
      
    sys.stderr.write("*** Reading cluster - barcode csv file ***\n\n")

    # Parse cluster file into dictionary
    clusters_dic = {}
    for line in open(barcode_annotations,'r'):
      if line.startswith("#"):
        continue
      line = line.rstrip().split(',')
      clusters_dic[line[1]] = line[0]

    # Iterate over clusters and generate list of paths for output files
    clusters = list(set(clusters_dic.values()))
    clusters_outfiles = {x: out_prefix + x + "_out.bam" for x in clusters}

    sys.stderr.write("*** Found following clusters in cluster - barcode file ***\n")
    sys.stderr.write("\n".join(clusters) + "\n\n")

    sys.stderr.write("*** Creating following output files ***\n")
    print("\n".join(clusters_outfiles.values()) + "\n\n")

    # Open bam file and save the header
    bamfile = pysam.AlignmentFile(bamFile, "rb")
    header  = bamfile.header

    # Iterate over bamfile
    # Write as bam

    with ExitStack() as stack:
        files = {fname: stack.enter_context(pysam.AlignmentFile(fname,'wb',header=header)) for fname in list(clusters_outfiles.values())}
        # Iterate over the bam file
        n = 0
        for line in bamfile:
          n+=1
          if n % 5000000 == 0:
            sys.stderr.write("*** {} lines processed ***\n".format(n))
          try:
            barcode = line.get_tag("CB")
          except KeyError:
            continue
          if barcode_prefix != "":
            barcode = barcode_prefix + "_" + barcode
          if barcode in clusters_dic:
            cluster = clusters_dic[barcode]
            files[clusters_outfiles[cluster]].write(line)


outdir_base = "data/"
bamBasePath = "/proj/tmp/mukund/atac/"


if not any(fname.endswith('_barcodes.csv') for fname in os.listdir(outdir_base)):
    barcodesToCsv(adata="data/210421_LOADED_merged_peak_matrix.h5ad", outdir="data")


for file_ in os.listdir(outdir_base):
    if file_.endswith("_barcodes.csv"):
        NGI_ID =  file_.split("_barcodes.csv")[0]
        OUTDIR = outdir_base + NGI_ID + "_bamFiles"
        if not os.path.isdir(OUTDIR):
            os.mkdir(OUTDIR)
            
        splitBam(bamFile=bamBasePath + NGI_ID +"/outs/possorted_bam.bam", 
        barcode_annotations=outdir_base+file_,
        out_prefix=OUTDIR+"/")
