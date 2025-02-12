#!/bin/sh
module load sambamba/0.6.7-gcc-4.8.5
module load boost/1.68.0-gcc-4.8.5
module load htslib/1.9-gcc-4.8.5
module load bcftools

echo "detect CNVs with $CNVdetector"
#bam is the first input parameter
bam=$1
CNVdetector="DELLY"
echo "detect CNVs with $CNVdetector"
#outDIR is the second input parameter
outDIR=$2

softpath=/BIGDATA2/scau_xlyuan_1/CJL/delly-main/src
excl=/BIGDATA2/scau_xlyuan_1/CJL/delly-main/excludeTemplates/human.hg19.excl.tsv
genome=/BIGDATA2/scau_xlyuan_1/xdt/WGS_WGBS_callingCNV/fasta/index/human/hg38.fa
sample=${bam/.sorted.bam/}
sample=${sample##*/}
echo $sample

${softpath}/delly call -x $excl \
-o $outDIR/${sample}.${CNVdetector}.bcf \
-g $genome $bam 
#bcftools code
bcftools view $outDIR/${sample}.${CNVdetector}.bcf > $outDIR/${sample}.${CNVdetector}.vcf
