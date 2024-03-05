

SAMPLE=$1 ## placeholder

#tar -xvf P27454_1003.tar.gz

echo "TAR extract complete"

#cd P27454_1003

#zcat P27454_1003.PT.pairs.gz | tail -n +401 > P27454_1003.PT.pairs

#tail -n +398 P27454_REdeep_hOPC.mapped.pairs > P27454_REdeep_hOPC.mapped.PT.pairs

echo "Valid pairs extracted! Creating .hic file..."

java -Xmx128g -Djava.awt.headless=true -jar ./juicer_tools.2.20.00.ac.jar pre --threads 64 P27454_REdeep_hOPC.mapped.PT.pairs P27454_REdeep_hOPC.mapped.pairs.contact-map.hic /data/proj/GCB_MK/micro-C/hg38.valid.chrom.sizes

echo "Created .hic file!"


#cd ..

#/home/mukund/miniconda3/envs/pycistopic/bin/python convertHic2Cool.py P27454_1003 5000

#echo "Created 5kb .cool file!"


## Run on merged hOPC replicates (1004, 1005)

#sort -k2,2V -k4,4V -k3,3n P27454_hOPC_merged.PT.pairs | \
#sed '1h;1d;$!H;$!d;G' > \
#P27454_hOPC_merged_sorted.PT.pairs

#echo "Sorted!"

#java -Xmx128g -Djava.awt.headless=true -jar juicer_tools.2.20.00.ac.jar pre --threads 64 P27454_hOPC_merged_sorted.PT.pairs P27454_hOPC_merged.PT.pairs.contact-map.hic /data/proj/GCB_MK/micro-C/hg38.valid.chrom.sizes
