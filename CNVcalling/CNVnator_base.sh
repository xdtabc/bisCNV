#!/bin/sh
set -e
module unload openssl
module load CNVnator
module load gcc/4.8.5
#module load root/6.12.04-gcc-4.8.5
##================##Predicting CNV regions##===============================##

bam=$1
outDIR=$2

bam_name=$(echo ${bam##*/})
ID=${bam_name/.bam/}
echo ${ID}

bam_name=`echo ${bam##*/}`
#echo ${bam_name}
echo 
rootname=${ID}.root

cnvnator -root ${rootname} -chrom \
chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY \
-tree $bam 
echo 'cnvnator -root [] -tree []' 

cnvnator -root ${rootname} -chrom \
chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY \
-his 100 -d /BIGDATA2/scau_xlyuan_1/xdt/WGS_WGBS_callingCNV/fasta/index/human/cnvnator
#用更名后的fa文件不行，
#参考基因组目录不需要作为参数
echo 'cnvnator -root [] -his 100 -d []' 
#PS：fasta's filename should be chr1.fa, chr2.fa, ... 
cnvnator -root ${rootname} -stat 100 
echo 'cnvnator -root [] -stat 100' 
cnvnator -root ${rootname} -partition 100 
echo 'cnvnator -root [] -partition 100' 
cnvnator -root ${rootname} -chrom \
chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY \
-call 100 > ${ID}.cnv
echo 'cnvnator -root [] -call 100 > []' 
export LC_ALL=C 
#先输入，不然perl显示LC_ALL=(unset) 
chmod a+x ${ID}.cnv
#code n
cnvnator2VCF.pl ${ID}.cnv > ${ID}.vcf 
#改名，加上cnvnator
mv ${ID}.vcf ${outDIR}/${ID}.cnvnator.vcf
	

