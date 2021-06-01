import os
import MOODS.parsers
import MOODS.tools
import MOODS.scan
from motif_utils import *


os.chdir("")

for chunk in chunk_bins:

    bin_GCdict = pickle.load(open(chunk["GCdict"], "rb")) ## 
    Peak = namedtuple('Peak', ['chrom', 'start', 'end', 'strand'])
    peak_iter = peak_reader(peak_name_list, select = [bin_GCdict["peaks"]])
    bg = bin_GCdict["counter"]
    for i in peak_iter:
        
    