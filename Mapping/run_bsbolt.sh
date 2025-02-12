#!/bin/sh

#bismark
#fastq1=$1
#fastq2=$2
#mappers=(bsbolt bwameth bismarkbt2 bsmap walt)
mapper="bsbolt"
echo "mapper: ${mapper}"
##============================================================##\
##============================================================##\
#运行4个bam 
#sample1-2-3-4 60x
samples=(sample1 sample2 sample3 sample4)
for j in {0..3}
do
    echo "$j"
    sample="${samples[j]}"
    #read_depths=(10x 15x 20x 30x 60x)
    read_depth="60x"
    echo "sample: ${sample}  RD: ${read_depth}"
    fastqDIR="/BIGDATA2/scau_xlyuan_1/xdt/WGS_WGBS_callingCNV/sherman_test/fastq/${sample}_${read_depth}"
    fastq1=$fastqDIR/"${sample}_${read_depth}_wgbs_BisCNV_R1.fastq.gz"
    fastq2=$fastqDIR/"${sample}_${read_depth}_wgbs_BisCNV_R2.fastq.gz"
    ls $fastq1
    yhbatch -J xdt${mapper:0:3}${read_depth} ${mapper}_pipeline.sh $fastq1 $fastq2
done
