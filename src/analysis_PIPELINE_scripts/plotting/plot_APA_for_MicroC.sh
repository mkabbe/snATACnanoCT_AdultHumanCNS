#!/bash


## This script runs the Aggregate Peak Analysis (pileup analysis) on the B cells and hOPCs
## 

conda activate hichip ## On Rackham

PAD=100 ## larger value generate larger "windows" for the pileup.

## APA analysis 

for loopfile in loops_mustache_*5kb.tsv
do
	echo "********* Running Loop ${loopfile}" # Track loop
	loopbase=$(echo $loopfile | cut -d_ -f3) #Store loop name
	cut -f1-6 $loopfile | tail -n +2 > TMP_loop.bedpe
	
	for coolfile in P27454_*5000.cool 
	do
		echo "****** Generating pileup on ${coolfile}"
		coolbase=$(echo $coolfile | cut -d. -f1 | cut -d_ -f2) # Store matrix name
		ofile=P27454_${coolbase}_mtx_v_${loopbase}_loop.txt

		coolpup.py $coolfile \
			TMP_loop.bedpe \
				--pad $PAD
				--unbalanced \
					--n_proc 24 \
						--outname $ofile
		
		echo "*** Completed ${ofile}"
	done
done



## Plot APA results

pileup_txt_files=$(ls P27454*_mtx_v_*loop.txt | tr "\n" " ")

plotpup.py $pileup_txt_files \
	--col_names BC_1.hic,BC_2.hic,BC_3.hic,hOPC.hic \
		--row_names BC_1_loop,BC_2_loop,BC_3_loop,hOPC_loop \
			--n_cols 4 \
				--vmin 1 --vmax 15 \
					--cmap coolwarm \
						--output P27454_Bcell_vs_hOPC_APA.pdf



