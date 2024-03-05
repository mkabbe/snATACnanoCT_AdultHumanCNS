#!/bin/python

import pickle 

## maps sampleID after cellranger-atac aggr to the NGI ID. 
## load the pickle file when populating anndata.obs 

sample_dict={}
ngi_id = "ref/NGI_id.txt"

with open(ngi_id) as f:
    for i,line in enumerate(f):
        sample_dict[str(i+1)] = line.strip()

with open("_sampleID.pickle", "wb") as handle:
    pickle.dump(sample_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)


