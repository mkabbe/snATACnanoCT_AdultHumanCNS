import gffutils
import pybedtools 
from pybedtools.featurefuncs import TSS
from gffutils.helpers import asinterval

#gffutils.create_db(data="ref/gencode.v35.annotation.gtf", dbfn="ref/gencode.v35.annotation.gtf.db")

db = gffutils.FeatureDB("ref/gencode.v35.annotation.gtf.db")

def tss_generator():
    """
        Generator function to yield TSS of each annotated transcript
        """
        for transcript in db.features_of_type('transcript'):
            yield TSS(asinterval(transcript), upstream=1, downstream=0)
    
# A BedTool made out of a generator, and saved to file.
tsses = pybedtools.BedTool(tss_generator()).saveas('ref/gencode_hg38_TSSes.gtf')

tsses_50bp = tsses.slop(b=25, genome='hg38', output='ref/tsses-50.gtf')

