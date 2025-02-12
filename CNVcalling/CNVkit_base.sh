#!/bin/sh
module unload anaconda3
module load anaconda3/2020.07
source activate cnvkit
#the env 
fasta="/BIGDATA2/scau_xlyuan_1/xdt/WGS_WGBS_callingCNV/fasta/index/human/hg38.fa"
bam=$1
CNVdetector="CNVkit"
echo "detect CNVs with $CNVdetector"
#outDIR is the second input parameter
outDIR=$2

dataDIR="/BIGDATA2/scau_xlyuan_1/xdt/WGS_WGBS_callingCNV/CNV_detect/script/cnv_code/CNVkit/data"
##===========##  target & antitarget
#bed用E:\wgs_wgbs_file\script\coverageMaster_mk_bed.py的代码创建
#输入bed准备好
name="hg38"
bed="${name}.area.bed"
refFlat="hg38.refFlat.txt"
#enter the directory of bed
bedDIR="/BIGDATA2/scau_xlyuan_1/xdt/WGS_WGBS_callingCNV/CNV_detect/script/cnv_code/CNVkit"
cd $outDIR
#每个物种是定值，不需生成第二次
#cnvkit.py target $bedDIR/$bed --annotate $bedDIR/${refFlat} -o ${bedDIR}/${name}.target.bed
#cnvkit.py antitarget ${name}.target.bed -g ${dataDIR}/access-5k-mappable.${name}.bed -o ${bedDIR}/${name}.antitarget.bed

### 下载软件
#wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/liftOver

##===========## 每个物种都要转换？ 
### 下载相关文件（hg38和hg19坐标转换）
#wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz 
#wget http://hgdownload.cse.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHg19.over.chain.gz

### 执行转换
#./liftOver ~/xdt/programs/cnvkit/data/access-5k-mappable.hg19.bed hg19ToHg38.over.chain.gz access-5k-mappable.hg38.bed unMapped.txt

#============##  mk cnn
ID=${bam##*/}
ID=${ID%%.*}
echo ${ID}
cnvkit.py coverage ${bam} ${bedDIR}/${name}.target.bed -o ${ID}.${name}.targetcoverage.cnn
cnvkit.py coverage ${bam} ${bedDIR}/${name}.antitarget.bed -o ${ID}.${name}.antitargetcoverage.cnn

#============##  reference
#定值，不需生成第二次
#cnvkit.py reference -o ${name}.FlatReference.cnn -f ${fasta} -t ${name}.target.bed -a ${name}.antitarget.bed

#===========##  fix
#sample_1_10x.chr1.antitargetcoverage.cnn为空文件
cnvkit.py fix \
${ID}.${name}.targetcoverage.cnn \
${ID}.${name}.antitargetcoverage.cnn \
${bedDIR}/${name}.FlatReference.cnn \
-o ${ID}.${name}.cnr

#===========##  segment
cnvkit.py segment ${ID}.${name}.cnr -o ${ID}.${name}.cns

#===========##  call
cnvkit.py call ${ID}.${name}.cns -o ${ID}.${name}.call.cns

##===========##  export
mkdir -p vcf 
cnvkit.py export vcf ${ID}.${name}.call.cns > ./vcf/${ID}.${name}.cnvkit.vcf