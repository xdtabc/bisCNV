#!/bin/sh
set -e 
module unload python
module unload openssl
module load samtools/1.11-gcc-4.8.5
module load bwa/0.7.17-gcc-4.8.5
module load gcc/4.8.5
module load gsl/2.1-gcc-4.8.5
module load hisat2/2.1.0-gcc-4.8.5
module load pigz

mapper="walt"
species="human"
indexDir="/BIGDATA2/scau_xlyuan_1/xdt/WGS_WGBS_callingCNV/fasta/index"/${species}/${mapper}/hg38.dbindex
fastqDir="/BIGDATA2/scau_xlyuan_1/xdt/WGS_WGBS_callingCNV/sherman_test/fastq"
fastq1_gz=$1
fastq2_gz=$2
sample=${fastq1_gz/1.fastq.gz/}
sample=${sample##*/}

sampleID=${sample/_wgbs_BisCNV_R/}
echo "sample: $sample  
      sampleID: $sampleID"

#generate the bam name
sam="${sample}.${mapper}.sam"
bam="${sample}.${mapper}.bam"
sortbam="${sample}.${mapper}.sorted.bam"

#the code of mapper
echo "map reads to ref genome using "$mapper
outDir="/BIGDATA2/scau_xlyuan_1/xdt/WGS_WGBS_callingCNV/pipeline/${species}/${mapper}/bam/${sampleID}"
time_Dir="/BIGDATA2/scau_xlyuan_1/xdt/WGS_WGBS_callingCNV/pipeline/${species}/${mapper}/time"
mkdir -p ${outDir} ${time_Dir} 
#this is important!
###==========================================================================================##
echo -e "mem\tRSS\trealTime\tcpusysTime\tcpuuserTime" > ${time_Dir}/${bam/.sorted.bam/}.Bench.csv
fastq1=${fastq1_gz/.gz/}
fastq2=${fastq2_gz/.gz/}
#the fastq must be in the format of txt
gzip -dc ${fastq1_gz} > ${fastq1}
gzip -dc ${fastq2_gz} > ${fastq2}
/usr/bin/time -f "%K\t%M\t%E\t%S\t%U" -o ${time_Dir}/${bam/.sorted.bam/}.Bench.csv -a \
/BIGDATA2/scau_xlyuan_1/gwt/software/walt-master/bin/walt \
-i ${indexDir} \
-t 8 \
-sam \
-1 ${fastq1} \
-2 ${fastq2} \
-o ${outDir}/${sam} && echo 'mapping is done' && rm ${fastq1} ${fastq2}

cd ${outDir}
#the code of samtools
samtools view -b -S ${sam} |\
samtools sort -@ 24 -o ${sortbam} && echo 'done' && rm ${sam}
samtools index -@ 24 ${sortbam} 
samtools flagstat ${sortbam} > ${sortbam}.flagstat
echo "mapping by $mapper is done"



##===================##use CNVdectors
CNVdetector_Names=(cn_mops CNVnator DELLY Pindel GASV CNVkit)

for ((i=0;i<${#CNVdetector_Names[@]};i++))
#for i in {6..6}
do 
    echo $i
    CNVdetector=${CNVdetector_Names[i]}
    echo ${CNVdetector}" is running"
    #µ÷ÓÃCNVdetector_base.sh½Å±¾
    echo "detect CNVs with $CNVdetector"
    #the base scripts only have two parameters 
    #bam is the first input parameter
    bam=$sortbam
    #outDIR is the second input parameter
    outDIR="/BIGDATA2/scau_xlyuan_1/xdt/WGS_WGBS_callingCNV/pipeline/${species}/${mapper}/${CNVdetector}/out"
    timeDIR="/BIGDATA2/scau_xlyuan_1/xdt/WGS_WGBS_callingCNV/pipeline/${species}/${mapper}/${CNVdetector}/time"
    codeDIR="/BIGDATA2/scau_xlyuan_1/xdt/WGS_WGBS_callingCNV/CNV_detect/script/cnv_code/${CNVdetector}"
    mkdir -p $outDIR $timeDIR
    echo "bamid: "${bam/.sorted.bam/}
    ls ${outDir}/$bam
    echo "check bam file"
    cd ${outDIR}/
    echo -e "mem\tRSS\trealTime\tcpusysTime\tcpuuserTime" >	${timeDIR}/${bam/.sorted.bam/}.${CNVdetector}.Bench.csv
    /usr/bin/time -f "%K\t%M\t%E\t%S\t%U" -o ${timeDIR}/${bam/.sorted.bam/}.${CNVdetector}.Bench.csv -a \
    bash ${codeDIR}/${CNVdetector}_base.sh ${outDir}/$bam $outDIR
done 