#!/bin/sh

#fastq1=$1
#fastq2=$2
#20x: sample3, sample4
mappers=(bsbolt bwameth bismarkbt2 bsmap walt)
for i in {4..4}
do
    mapper="${mappers[i]}"
    echo "mapper: ${mapper}"
    samples=(sample1 sample2 sample3 sample4)
    for j in {0..0}
    do
        echo "$j"
        sample="${samples[j]}"
        read_depths=(10x 15x 20x 30x 60x)
        for k in {4..4}
        do
            read_depth=${read_depths[k]}
            fastqDIR="/BIGDATA2/scau_xlyuan_1/xdt/WGS_WGBS_callingCNV/sherman_test/fastq/${sample}_${read_depth}"
            fastq1=$fastqDIR/"${sample}_${read_depth}_wgbs_BisCNV_R1.fastq.gz"
            fastq2=$fastqDIR/"${sample}_${read_depth}_wgbs_BisCNV_R2.fastq.gz"
            yhbatch -J xdt${mapper:0:3}${read_depth} ${mapper}_pipeline.sh $fastq1 $fastq2
        done
    done
done

#30x
#mappers=(bsbolt bwameth bismarkbt2 bsmap walt)
#for i in {4..4}
#do
#    mapper="${mappers[i]}"
#    echo "mapper: ${mapper}"
#    samples=(sample1 sample2 sample3 sample4)
#    for j in {0..3}
#    do
#        echo "$j"
#        sample="${samples[j]}"
#        read_depths=(10x 15x 20x 30x 60x)
#        for k in {3..3}
#        do
#            read_depth=${read_depths[k]}
#            fastqDIR="/BIGDATA2/scau_xlyuan_1/xdt/WGS_WGBS_callingCNV/sherman_test/fastq/${sample}_${read_depth}"
#            fastq1=$fastqDIR/"${sample}_${read_depth}_wgbs_BisCNV_R1.fastq.gz"
#            fastq2=$fastqDIR/"${sample}_${read_depth}_wgbs_BisCNV_R2.fastq.gz"
#            yhbatch -J xdt${mapper:0:3}${read_depth} ${mapper}_pipeline.sh $fastq1 $fastq2
#        done
#    done
#done

#60x fastq需要补充
