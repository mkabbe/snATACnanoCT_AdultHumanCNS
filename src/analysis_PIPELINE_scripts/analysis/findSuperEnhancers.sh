#conda activate rose_py2
cd ../rose_SE/rose/

base_dir="/data/proj/GCB_MK/scCT/nanoCT_EAE/fragments/P29054/ctype/H3K27ac/"
for ctype in AST CXEX CXINH #ENDO
do

    cp ${base_dir}/NARROWPEAKS/sorted_${ctype}_H3K27ac_fragments.bed_peaks.narrowPeak ${base_dir}/NARROWPEAKS/sorted_${ctype}_H3K27ac_peaks.bed

    python ROSE_main.py \
    -g HG38 \
    -i ${base_dir}/NARROWPEAKS/sorted_${ctype}_H3K27ac_peaks.bed \
    -r ${base_dir}/BW_FILES/${ctype}_H3K27ac_fragments.bed.sorted.bam \
    -o ${ctype}_roseSE \
    -t 3000

done

## call region-specific super-enhancers
# base_dir="/data/proj/GCB_MK/scCT/nanoCT_EAE/fragments/P29054/OLG_tissue/H3K27ac/"
# for ctype in OLG_CSC OLG_BA4
# do

#     cp ${base_dir}/NARROWPEAKS/sorted_${ctype}_H3K27ac_fragments.bed_peaks.narrowPeak ${base_dir}/NARROWPEAKS/sorted_${ctype}_H3K27ac_peaks.bed
    
#     python ROSE_main.py \
#     -g HG38 \
#     -i ${base_dir}/NARROWPEAKS/sorted_${ctype}_H3K27ac_peaks.bed \
#     -r ${base_dir}/BW_FILES/${ctype}_H3K27ac_fragments.bed.sorted.bam \
#     -o ${ctype}_roseSE \
#     -t 3000

# done