find_names<-function(patname, path){
  #需要输入后缀名
  #返回符合后缀名的文件
  outnames=c()
  #当前目录下查找文件
  #path=getwd()
  for (name in dir(path)){
    #name='sample1.cn_mops.bwa.bed'
    a=strsplit(name, split = '[.]')[[1]]
    #break}
    if (patname %in% a){
      #如果包含指定子字符串，返回T；否则返回F
      #添加到一个向量中
      outnames<-c(outnames, name)
    }else{
      next
    }
  }
  return(outnames)
}

#==============================================================================#
#recall
#0324
args <- commandArgs(trailingOnly = T)
print(args[1])
indir=args[1]
setwd(indir)
#dir()
#setwd(outDIR)
#cn_mops
patname='cn_mops'
filenames=find_names(patname, getwd())
print(34)
print(filenames)
#开始统计
firstDF=T
for (bedname in filenames){
  bed=read.delim(bedname,header = F)
  head(bed)
  group<-gsub('.overlap', '', bedname)
  #
  colnames(bed)=c('chr2','start2','end2','cnvlen2','cnvtype2','cnvname',
                  'chr1','start1','end1','cnvlen1','cnvtype1','overlaplen')
  head(bed)
  #统计每个bed中有多少simuCNV
  simuCNVnum=length(unique(paste(bed$chr1,bed$start1)))
  table(bed$cnvtype1)
  sort(unique(bed$cnvtype1))
  sort(unique(bed$cnvtype2))
  #重新命名
  #现在进行筛选
  #给bed$cnvtype1重新赋值便于匹配
  #
  unique(bed$CNVtype)
  head(bed)
  SUM1<-data.frame(matrix(0,2,5))
  #修改列名
  colnames(SUM1)<-c('group','cnvtype','positiveNUM','CNVNUM','positiveRATE')
  #对应关系
  ###=================================================###
  ##CNVER的CNVTYPE
  bed$CNVtype1=0
  sort(unique(bed$cnvtype1))
  #del
  NAMEs=c('DEL')
  for (i in 1:length(NAMEs)){
    NAME=NAMEs[i]
    #赋值为del
    bed$CNVtype1[grep(NAME,bed$cnvtype1)]='del'
  }
  #dup
  NAMEs=c('DUP')
  for (i in 1:length(NAMEs)){
    NAME=NAMEs[i]
    #赋值为del
    bed$CNVtype1[grep(NAME,bed$cnvtype1)]='dup'
  }
  ###=================================================###
  ##REFCNV的CNVTYPE
  bed$CNVtype2=0
  NAMEs=c('del')  #del不变
  for (i in 1:length(NAMEs)){
    NAME=NAMEs[i]
    #赋值为del
    bed$CNVtype2[grep(NAME,bed$cnvtype2)]='del'
  }
  #dup
  NAMEs=c("duporigin","insorigin")
  for (i in 1:length(NAMEs)){
    NAME=NAMEs[i]
    #赋值为del
    bed$CNVtype2[grep(NAME,bed$cnvtype2)]='dup'
  }
  ###=================================================###
  head(bed)
  cnvtypes=c('del','dup')
  #生成相同的名称了
  for (i in 1:2){
    #i=1
    #选出所有相同类型的softCNV
    df1=bed[bed$CNVtype1==cnvtypes[i],]
    #去重了，recall手动计算分母，每个物种不同!!
    #sample(1-4).simu.bed
    dirname="simu_human"
    adir=paste("/BIGDATA2/scau_xlyuan_1/xdt/WGS_WGBS_callingCNV/pipeline/CNV_result/refCNV/",dirname,'/',sep="")
    #判别当前df1属于哪个BED
    #取与bedname一致的sample.simu.bed
    simubed_name=paste(adir,unlist(strsplit(bedname, split = '_'))[1], ".simu.bed", sep="")
    print(simubed_name)
    simubed=read.table(simubed_name)
    if (i==1){
      detectCNVNUM=nrow(simubed[simubed$V5=="del",])
    }else{
      detectCNVNUM=nrow(simubed[simubed$V5=="duporigin"|simubed$V5=="insorigin",])
    }
    #df1中全是该CNV类型
    b=df1[df1$CNVtype2==cnvtypes[i],]
    b
    #去重
    positiveNUM=length(unique(paste(b$chr1,b$start1,b$end1)))
    #算比率
    positiveRATE=positiveNUM/detectCNVNUM
    #写入表格
    SUM1[i,1]<-group
    SUM1[i,2]<-cnvtypes[i]
    SUM1[i,3]<-positiveNUM
    SUM1[i,4]<-detectCNVNUM
    SUM1[i,5]<-positiveRATE
  }
  if (firstDF){
    rateSUM<-SUM1
    #关闭开关
    firstDF<-F
  }else{
    rateSUM<-rbind(rateSUM, SUM1)
  }
}
print(130)
rateSUM
rateSUM$group=gsub(".overlap|.sorted|.bam","",rateSUM$group)
rateSUM<-rateSUM[order(rateSUM$group),]
print(args[2])
#outname包括路径
outname=args[2]
write.csv(rateSUM, outname, row.names = F)
