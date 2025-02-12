#!/bin/sh
module unload anaconda3
module load anaconda3/2020.07
source activate py27
module load samtools
fa="/BIGDATA2/scau_xlyuan_1/xdt/WGS_WGBS_callingCNV/fasta/index/human/hg38.fa"
#bam is the first input parameter
bam=$1
CNVdetector="Pindel"
echo "detect CNVs with $CNVdetector"
#outDIR is the second input parameter
outDIR=$2

echo "outDIR: $outDIR"
DATE=$(date +%Y%m%d)
bam_name=${bam/.sorted.bam/}
echo "17"
bam_name=${bam_name##*/}
echo "bam_name: $bam_name"
mkdir -p ${outDIR}/${bam_name%%.*}
cd ${outDIR}/${bam_name%%.*}
echo ${outDIR}/${bam_name%%.*}

echo "23"
configFILE=${outDIR}/${bam_name%%.*}/config.txt
#insert size: 250
echo ${bam}' 250 '${bam_name%%.*} > ${configFILE}
pindel -T 20 -f ${fa} -i ${configFILE} -o pindel.out
outnames=(find ${outDIR} -name 'pindel.out_*')
echo "outnames: "${outnames[@]}
for outfile in ${outnames[@]}
do
	pindel2vcf \
	-r ${fa} \
	-R hg38 \
	-p ${outfile} \
	-d ${DATE} \  # 
	-v ${outfile}.vcf 
done
echo "39"
cat pindel.out_*.vcf > ${bam_name%%.*}.vcf 