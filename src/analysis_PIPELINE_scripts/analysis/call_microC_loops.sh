
#cooler balance P27454_hOPC_REdeep/P27454_REdeep_hOPC.mapped.pairs.contact-map.25000.cool
#cooler balance P27454_hOPC_REdeep/P27454_REdeep_hOPC.mapped.pairs.contact-map.100000.cool
#cooler balance obs_exp_hOPC_chr17_5000.cool --cis-only 

#conda activate hicexplorer

#cd P27454_1001
#mustache -f P27454_1001.PT.pairs.contact-map.5000.cool -r 5kb -pt 0.001 -o loops_mustache_1001_5kb.tsv


#cd ../P27454_1002

#cooler balance P27454_1002.PT.pairs.contact-map.5000.cool
#mustache -f P27454_1002.PT.pairs.contact-map.5000.cool -r 5kb -pt 0.001 -o loops_mustache_1002_5kb.tsv


#cd ../P27454_1003

#cooler balance P27454_1003.PT.pairs.contact-map.5000.cool
#mustache -f P27454_1003.PT.pairs.contact-map.5000.cool -r 5kb -pt 0.001 -o loops_mustache_1003_5kb.tsv


cd ./P27454_hOPC_REdeep

#cooler balance P27454_REdeep_hOPC.mapped.pairs.contact-map.5000.cool
mustache -f P27454_REdeep_hOPC.mapped.pairs.contact-map.5000.cool -r 5kb -pt 0.01 -o loops_long_mustache_hOPC_5kb.tsv --distance 100000000

