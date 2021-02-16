#!/usr/bin/env python

import sys 

def binGenome(binsize, stepsize = 0, chromsize_file, outfile):
    '''
    binsize:  feature window size
    stepsize: overlap size for adjacent windows
    '''

    bin_index = 1
    binsize=int(binsize)
    stepsize=int(stepsize)
    with open(outfile,"w") as out:
        with open(chromsize_file, "r") as f_:
            if stepsize == 0:
                stepsize=binsize

            for line in f_:
                line = line.strip().split("\t")
                chr_name = line[0]
                chr_len = int(line[1])

                startpos = 0
                #open(outfile, "w")
                for endpos in range(binsize, chr_len + stepsize, stepsize):
                    if min(chr_len, endpos) > startpos:
                        out.write("{0}\t{1}\t{2}\t{3}\n".format(chr_name, startpos, min(endpos,chr_len),bin_index))
                        startpos += stepsize
                        bin_index += 1


binGenome(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])

