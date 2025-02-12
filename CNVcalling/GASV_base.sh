#!/bin/sh
#yhbatch -N 1 -p rhenv -J xdtgasv gasv_0311.sh
echo "scriptNAME: $0"
module load samtools/1.11-gcc-4.8.5
module load java

#bam pretreatment
bam=$1
outDIR=$2
#cd $outDIR
echo "11"
ls $bam $outDIR
java -Xms512m -Xmx2048m \
-jar /BIGDATA2/scau_xlyuan_1/CJL/gasv/bin/BAMToGASV.jar ${bam} #-MAPPING_QUALITY 30 -CUTOFF_LMINLMAX SD=50
#java -jar /BIGDATA2/scau_xlyuan_1/CJL/gasv/bin/BAMToGASV.jar ${bam}

echo "31"
java -jar /BIGDATA2/scau_xlyuan_1/CJL/gasv/bin/GASV.jar \
--batch ${bam}.gasv.in
ls -l ${bam}.gasv.in
echo "34"
echo "outname: ${bam}.gasv.in.clusters"
mv ${bam}.gasv.in.clusters $outDIR
