#!/bin/bash

### Plot metasignal for all HOX clusters for OLs and OPCs split by tissue

computeMatrix reference-point -R hox_regions/hoxa_promoters.bed -S BW_FILES/*bw -a 2000 -b 2000 -o hox_regions/hoxa_coverage.matrix
computeMatrix reference-point -R hox_regions/hoxb_promoters.bed -S BW_FILES/*bw -a 2000 -b 2000 -o hox_regions/hoxb_coverage.matrix
computeMatrix reference-point -R hox_regions/hoxc_promoters.bed -S BW_FILES/*bw -a 2000 -b 2000 -o hox_regions/hoxc_coverage.matrix
computeMatrix reference-point -R hox_regions/hoxd_promoters.bed -S BW_FILES/*bw -a 2000 -b 2000 -o hox_regions/hoxd_coverage.matrix



plotProfile -m hox_regions/hoxa_coverage.matrix -o hox_regions/plots/olg_hoxa_coverage.pdf -T HOXA_coverage --perGroup --plotHeight 20 --plotWidth 30
plotProfile -m hox_regions/hoxb_coverage.matrix -o hox_regions/plots/olg_hoxb_coverage.pdf -T HOXB_coverage --perGroup --plotHeight 20 --plotWidth 30
plotProfile -m hox_regions/hoxc_coverage.matrix -o hox_regions/plots/olg_hoxc_coverage.pdf -T HOXC_coverage --perGroup --plotHeight 20 --plotWidth 30
plotProfile -m hox_regions/hoxd_coverage.matrix -o hox_regions/plots/olg_hoxd_coverage.pdf -T HOXD_coverage --perGroup --plotHeight 20 --plotWidth 30