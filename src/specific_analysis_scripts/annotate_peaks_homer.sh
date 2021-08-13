#!/bin/bash

## HOMER annotation of differentially accessible peaks


annotatePeaks.pl sig_diff_peaks_OPC_CSCvsBA4.bed hg38 | sort -k1,1 -V > annot_sig_diff_peaks_OPC_CSCvsBA4.bed
echo -e "logfoldchange \tpval \tpval_adj" > tmp && cut -f6,7,8 sig_diff_peaks_OPC_CSCvsBA4.bed >> tmp
cut -f6,7,8 sig_diff_peaks_OPC_CSCvsBA4.bed >> tmp
paste annot_sig_diff_peaks_OCP_CSCvsBA4.bed tmp > annot && mv annot annot_sig_diff_peaks_OPC_CSCvsBA4.bed