#!/bin/sh
set -e 
module unload python
module unload openssl
module unload samtools
module load samtools/1.9-gcc-4.8.5
module load bwa/0.7.17-gcc-4.8.5
module load gcc/4.8.5
module load gsl/2.1-gcc-4.8.5
module load hisat2/2.1.0-gcc-4.8.5
module load pigz

mapper="bwameth"
species="human"
indexDir="/BIGDATA2/scau_xlyuan_1/xdt/WGS_WGBS_callingCNV/fasta/index"/${species}/${mapper}/hg38.fa
ls ${indexDir}
fastqDir="/BIGDATA2/scau_xlyuan_1/xdt/WGS_WGBS_callingCNV/sherman_test/fastq"
#sample="ref_test_BisCNV_R"
fastq1=$1
fastq2=$2
ls ${fastq1}
sample=${fastq1/1.fastq.gz/}
sample=${sample##*/}
echo "$sample"
sampleID=${sample/_wgbs_BisCNV_R/}

echo "sample: $sample  
      sampleID: $sampleID"
sam="${sample}.${mapper}.sam"
bam="${sample}.${mapper}.bam"
sortbam="${sample}.${mapper}.sorted.bam"

echo "map reads to ref genome using "$mapper

#the code of the mapper
outDir="/BIGDATA2/scau_xlyuan_1/xdt/WGS_WGBS_callingCNV/pipeline/${species}/${mapper}/bam/${sampleID}"
time_Dir="/BIGDATA2/scau_xlyuan_1/xdt/WGS_WGBS_callingCNV/pipeline/${species}/${mapper}/time"
mkdir -p ${outDir} ${time_Dir} 
echo -e "mem\tRSS\trealTime\tcpusysTime\tcpuuserTime" >	${time_Dir}/${bam/.sorted.bam/}.Bench.csv
/usr/bin/time -f "%K\t%M\t%E\t%S\t%U" -o ${time_Dir}/${bam/.sorted.bam/}.Bench.csv -a \
/BIGDATA2/scau_xlyuan_1/gwt/software/bwa-meth-master/bwameth.py \
--reference ${indexDir} \
${fastq1} \
${fastq2} \
-t 8 \
> ${outDir}/${sam}

cd ${outDir}
#the code of samtools
ls -l ${sam}
echo ${sam}

samtools view -b -S ${sam} -o ${bam} && rm ${sam}
samtools sort ${bam} -@ 24 -o ${sortbam} && echo 'done' && rm ${bam}
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
