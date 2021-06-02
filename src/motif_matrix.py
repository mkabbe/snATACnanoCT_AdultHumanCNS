import os
import MOODS.parsers
import MOODS.tools
import MOODS.scan
from motif_utils import Motifs


os.chdir("")

for chunk in chunk_bins:

    bin_GCdict = pickle.load(open(chunk["GCdict"], "rb")) ## 
    Peak = namedtuple('Peak', ['chrom', 'start', 'end', 'strand'])
    peak_iter = peak_reader(peak_name_list, select = [bin_GCdict["peaks"]])
    bg = bin_GCdict["counter"]
    for i in peak_iter:

## convert jaspar to PFM 
jaspar_path = "pfm_files/"
for jaspar_file in os.listdir(jaspar_path):
    if jaspar_file.endswith(".jaspar"):
        jpath = jaspar_path.rstrip("/") + "/" ## handles missing/trailing "/"
        Motifs.jaspar_to_pfm(f"{jpath}/{jaspar_file}")


peak_fa_list=[]
peak_name_list=[]
for header,seq in iter_fasta(PEAK_FASTA):
    peak_fa_list.append(seq) # store the FASTA sequences in a list
    peak_name_list.append(header) # store the FASTA headers in a list (chr:start-stop)
    
    
                    



'''
convert jaspar to pfm

for bin in GCbin:
        load the bin-Dict.
        get bg value
        get Peak list
        get FASTA sequences
        
        make matrices, thresholds (using bg value) from pfm
        
        scan_motif_from_bed()
        

'''




def greet(x):
    print(f"Greetings, {x}!")
    
class Test:
    def __init__(self, name, age):
        self.name = name
        self.age = age
    def getName(self, new):
        print(f'Name is {self.name}')
        print(f'New name is {new}')
        greet(self.name)
        greet(new)
    def getAge(self):
        print(f'Age is {self.age}')
    