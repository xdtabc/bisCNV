#!/bin/sh
#确认breakdancer-max能运行

#bam is the first input parameter
bam=$1
#outDIR is the second input parameter
outDIR=$2
ID=${bam##*/}
ID=${ID/.sorted.bam/}
echo "ID: "${ID}
echo "11"
echo "$outDIR"

#right now it's not in any conda env
/home/wenping/miniconda3/pkgs/breakdancer-1.4.5-0/lib/breakdancer-maxunstable/bam2cfg.pl \
${bamdir}/${bam} -g -h > ${outDIR}/${ID}.cfg

#the env is important!!
source activate xdt
breakdancer-max -q 10 ${outDIR}/${ID}.cfg > ${outDIR}/${ID}.ctx 
