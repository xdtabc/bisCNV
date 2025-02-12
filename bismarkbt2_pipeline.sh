#!/bin/sh
#conda activate bsbolt
export LC_ALL=C
module load bwa/0.7.17-gcc-4.8.5
module unload python
module load python/2.7.9-gcc-4.8.5
module unload samtools 
module unload openssl
module load samtools/1.11-gcc-4.8.5
module load java

mapper="bismarkbt2"
species="human"
indexDir="/BIGDATA2/scau_xlyuan_1/xdt/WGS_WGBS_callingCNV/fasta/index"/${species}/${mapper}
fastqDir="/BIGDATA2/scau_xlyuan_1/xdt/WGS_WGBS_callingCNV/sherman_test/fastq"
#sample="ref_test_BisCNV_R"
fastq1=$1
fastq2=$2
ls ${fastq1}

#fastq1="sample1_10x_wgbs_BisCNV_R1.fastq.gz"
sample=${fastq1/1.fastq.gz/}
sample=${sample##*/}
sampleID=${sample/_wgbs_BisCNV_R/}
echo "sample: $sample  
      sampleID: $sampleID"
#sample="sample1_10x_wgbs_BisCNV_R"
bam="${sample}1_bismark_bt2_pe.bam"
sortbam="${sample}.${mapper}.sorted.bam"

echo "map reads to ref genome using $mapper"
#the code of the mapper
outDir="/BIGDATA2/scau_xlyuan_1/xdt/WGS_WGBS_callingCNV/pipeline/${species}/${mapper}/bam/${sampleID}"
time_Dir="/BIGDATA2/scau_xlyuan_1/xdt/WGS_WGBS_callingCNV/pipeline/${species}/${mapper}/time"
tempDir="/BIGDATA2/scau_xlyuan_1/xdt/WGS_WGBS_callingCNV/pipeline/${species}/${mapper}/temp"
mkdir -p ${outDir} ${time_Dir} ${tempDIR}
echo "35"
echo -e "mem\tRSS\trealTime\tcpusysTime\tcpuuserTime" >	${time_Dir}/${bam/.sorted.bam/}.Bench.csv
/usr/bin/time -f "%K\t%M\t%E\t%S\t%U" -o  ${time_Dir}/${bam/.sorted.bam/}.Bench.csv -a \
time /BIGDATA2/scau_xlyuan_1/gwt/software/bismark/Bismark-0.22.3/bismark \
--path_to_bowtie2 /BIGDATA2/scau_xlyuan_1/gwt/software/bowtie2/bowtie2-2.3.5.1-linux-x86_64 \
--genome ${indexDir} \
-1 ${fastq1} \
-2 ${fastq2} \
--parallel 6 \
${tempDir} \
-o ${outDir} \
--ambiguous --ambig_bam

echo "48"
cd ${outDir}/
sampleID="${sample/_wgbs_BisCNV_R/}"
mkdir -p ${outDir}/${sampleID}_ambiguous

mv ${sample}1.fastq.gz_bismark_bt2_pe.ambig.bam ${sample}1_bismark_bt2_PE_report.txt ${outDir}/${sampleID}_ambiguous
rm ${sample}1.fastq.gz_ambiguous_reads_1.fq.gz ${sample}2.fastq.gz_ambiguous_reads_2.fq.gz
#排序
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
