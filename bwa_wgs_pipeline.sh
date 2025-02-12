#!/bin/sh
#conda activate bsbolt
set -e 
module unload openssl
module load samtools/1.11-gcc-4.8.5
module load java
module load bwa/0.7.17-gcc-4.8.5
module load anaconda3/2020.07
source activate BSBolt

echo "11"
#WGS
mapper="bwa"
species="human"
genomeID="hg38"
indexDir="/BIGDATA2/scau_xlyuan_1/xdt/WGS_WGBS_callingCNV/fasta/index"/${species}/${mapper}/${genomeID}.fa
ls ${indexDir}

fastq1=$1
fastq2=$2
sample=${fastq1/1.fastq.gz/}
sample=${sample##*/}
echo "$sample"
sampleID=${sample/_wgbs_BisCNV_R/}
echo "sample: $sample  
      sampleID: $sampleID"
sam="${sample}.${mapper}.sam"
bam="${sample}.bam"
sortbam="${sample}.${mapper}.sorted.bam"

echo "map reads to ref genome using $mapper"
#the code of the mapper
outDir="/BIGDATA2/scau_xlyuan_1/xdt/WGS_WGBS_callingCNV/pipeline/${species}/${mapper}/bam"
time_Dir="/BIGDATA2/scau_xlyuan_1/xdt/WGS_WGBS_callingCNV/pipeline/${species}/${mapper}/time"
mkdir -p ${outDir} ${time_Dir} 
#this is important!
###==========================================================================================##
echo -e "mem\tRSS\trealTime\tcpusysTime\tcpuuserTime" > ${time_Dir}/${bam/.sorted.bam/}.Bench.csv
/usr/bin/time -f "%K\t%M\t%E\t%S\t%U" -o ${time_Dir}/${bam/.sorted.bam/}.Bench.csv -a \

bwa mem -t 24 -M -R '@RG\tID:group_n\tLB:library_n\tPL:illumina\tPU:unit1\tSM:sample_n' \
${indexDir} ${fastq1} ${fastq2} |samtools view -b -S -o ${bam} 
echo "${outDir}"

cd ${outDir}/
ls $bam
samtools sort -@ 24 ${bam} -o $sortbam && rm ${bam}
samtools index $sortbam
samtools flagstat $sortbam > $sortbam.flagstat
echo "mapping by $mapper is done"


##===================##use CNVdectors
CNVdetector_Names=(cn_mops CNVnator DELLY Pindel GASV CNVkit)

for ((i=0;i<${#CNVdetector_Names[@]};i++))
#for i in {6..6}
do 
    echo $i
    CNVdetector=${CNVdetector_Names[i]}
    echo ${CNVdetector}" is running"
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

