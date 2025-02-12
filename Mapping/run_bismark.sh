#!/bin/sh

#bismark
#fastq1=$1
#fastq2=$2
#mappers=(bsbolt bwameth bismarkbt2 bsmap walt)
mapper="bismarkbt2"
echo "mapper: ${mapper}"
#samples=(sample1 sample2 sample3 sample4)
#60x sample1 sample3
#60x sample4
#samples=(sample3)
#for sample in ${samples[@]}
#do
#    read_depths=(10x 15x 20x 30x 60x)
#    read_depth="60x"
#    echo "sample: ${sample} read_depth: ${read_depth}"
#    fastqDIR="/BIGDATA2/scau_xlyuan_1/xdt/WGS_WGBS_callingCNV/sherman_test/fastq/${sample}_${read_depth}"
#    fastq1=$fastqDIR/"${sample}_${read_depth}_wgbs_BisCNV_R1.fastq.gz"
#    fastq2=$fastqDIR/"${sample}_${read_depth}_wgbs_BisCNV_R2.fastq.gz"
#    ls $fastq1
#    yhbatch -J xdt${mapper:0:3}${read_depth} ${mapper}_pipeline.sh $fastq1 $fastq2
#done

species="human"
dataTYPE="WGBS"
echo "species: ${species}   dataTYPE: ${dataTYPE}"
samples=(sample1 sample2 sample3)
for sample in ${samples[@]}
do
    echo "sample: ${sample}"
    fastqDIR="/BIGDATA2/scau_xlyuan_1/xdt/WGS_WGBS_callingCNV/realDATA/fastq/${species}/${dataTYPE}/${sample}"
    fastqs=(`ls $fastqDIR/*1.fastq.gz`)
    for fastq in ${fastqs[@]}
    do
        fqID=${fastq##*/}
        fqID=${fqID/1.fastq.gz/}
        echo $fqID
        fastq1=$fastqDIR/"${fqID}1.fastq.gz"
        fastq2=$fastqDIR/"${fqID}2.fastq.gz"
        ls $fastq1
        yhbatch -J xdt${mapper:0:3}${read_depth} ${mapper}_pipeline.sh $fastq1 $fastq2
    done
done
