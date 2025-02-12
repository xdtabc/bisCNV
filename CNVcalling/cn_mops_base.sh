#!/bin/sh
#conda activate xdt
module load R/3.6.3-gcc-4.8.5
ID=${sample}
echo ${ID}
CNVdetector="cn_mops"
echo "detect CNVs with $CNVdetector"
#bam is the first input parameter
bam=$1
#outDIR is the second input parameter
outDIR=$2
codeDIR="/BIGDATA2/scau_xlyuan_1/xdt/WGS_WGBS_callingCNV/CNV_detect/script/cnv_code/${CNVdetector}"

#soft code
Rscript ${codeDIR}/cnmops.R ${bam} ${outDIR}
