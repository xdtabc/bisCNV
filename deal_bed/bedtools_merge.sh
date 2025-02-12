#!/bin/sh

module load bedtools2/2.26.0-gcc-4.8.5


baseDIR=/BIGDATA2/scau_xlyuan_1/xdt/WGS_WGBS_callingCNV/CNV_detect/overlap/filter_bed
bedDIR=${baseDIR}/simu_1
outDIR=${baseDIR}/simu_2
#bedDIR=${baseDIR}/real_1
#outDIR=${baseDIR}/real_2
mkdir -p ${outDIR}


cd ${bedDIR}

#beds=(`find ${bedDIR} -name "*.bed"|sort`)
#
#for bed in ${beds[@]}
#do
#echo ${bed}
#mergebed=${bed/.filter.1.bed/}.merge.bed
#bedtools merge -i ${bed} -c 5 -o distinct > ${mergebed}
#mv ${mergebed} ${outDIR}
#
#done


cd /BIGDATA2/scau_xlyuan_1/xdt/WGS_WGBS_callingCNV/CNV_detect/overlap/filter_bed/all_ref_bed
bed=real_CNV.bed.dupNdel.txt
mergebed=real_CNV.dupNdel.merge.bed
bedtools merge -i ${bed} -c 6 -o distinct > ${mergebed}