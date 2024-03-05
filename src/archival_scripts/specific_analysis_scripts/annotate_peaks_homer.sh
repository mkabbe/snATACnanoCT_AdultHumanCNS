#!/bin/bash

## HOMER annotation of differentially accessible peaks
## run on Monod
## CL argument is the BED file ($1)

annotatePeaks.pl $1 hg38 | sort -k1,1 -V > annot_$1
echo -e "logfoldchange \tpval \tpval_adj" > tmp && cut -f6,7,8 $1 >> tmp
paste annot_$1 tmp > annot && mv annot annot_$1

