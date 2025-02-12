#!/bin/sh

module load bedtools2/2.26.0-gcc-4.8.5


baseDIR=/BIGDATA2/scau_xlyuan_1/xdt/WGS_WGBS_callingCNV/CNV_detect/overlap/filter_bed
/BIGDATA2/scau_xlyuan_1/xdt/WGS_WGBS_callingCNV/CNV_detect/overlap/filter_bed/
#bedDIR=${baseDIR}/simu_2
bedDIR=${baseDIR}/real_2
#outDIR=${baseDIR}/simu2true_overlap
outDIR=${baseDIR}/real2true_overlap
mkdir -p ${outDIR}


cd ${bedDIR}



refDIR=/BIGDATA2/scau_xlyuan_1/xdt/WGS_WGBS_callingCNV/CNV_detect/overlap/filter_bed/all_ref_bed

#afile=${refDIR}/simu.all.merge.bed
afile=${refDIR}/real_CNV.dupNdel.1.merge.bed

bfiles=(`find ${bedDIR} -name '*.bed'|sort`)

for bfile in ${bfiles[@]}
do
	echo ${bfile}

	bfileID=`echo ${bfile}|cut -d / -f 10`
	bfileID=${bfileID/.filter.bed.merge.bed/}
	echo ${outDIR}/${bfileID}
	bedtools intersect -a ${afile} -b ${bfile} -c > ${outDIR}/${bfileID}.overlap
done

bfile=/BIGDATA2/scau_xlyuan_1/xdt/WGS_WGBS_callingCNV/CNV_detect/overlap/filter_bed/real_2/realDATA_WGS.cnvnator.merge.bed
	bfileID=`echo ${bfile}|cut -d / -f 10`
	bfileID=${bfileID/.filter.bed.merge.bed/}
bedtools intersect -a ${afile} -b ${bfile} -c > ${outDIR}/${bfileID}.overlap

#realign
#precision
module load bedtools2/2.26.0-gcc-4.8.5
outDIR1=/BIGDATA2/scau_xlyuan_1/xdt/WGS_WGBS_callingCNV/WGBS_REalign/cnvout
afile=/BIGDATA2/scau_xlyuan_1/xdt/WGS_WGBS_callingCNV/WGBS_REalign/cnvout/sherman_sample_3_10x_WGBS_simulated_bsmap.sorted_changed.1.sorted.1.cnvnator.bed
bfile=/BIGDATA2/scau_xlyuan_1/xdt/WGS_WGBS_callingCNV/CNV_detect/overlap/simu_ref_cnv_bed/sample3.simu.bed
bedtools intersect -a ${afile} -b ${bfile} -wao > ${outDIR1}/${bfileID}.precision.overlap
#recall
afile=/BIGDATA2/scau_xlyuan_1/xdt/WGS_WGBS_callingCNV/CNV_detect/overlap/simu_ref_cnv_bed/sample3.simu.bed
bfile=/BIGDATA2/scau_xlyuan_1/xdt/WGS_WGBS_callingCNV/WGBS_REalign/cnvout/sherman_sample_3_10x_WGBS_simulated_bsmap.sorted_changed.1.sorted.1.cnvnator.bed
bedtools intersect -a ${afile} -b ${bfile} -wao > ${outDIR1}/${bfileID}.recall.overlap

