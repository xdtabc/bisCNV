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
#precision
#0324
args <- commandArgs(trailingOnly = T)
print(args[1])
indir=args[1]
setwd(indir)
#dir()
#setwd(outDIR)
patname='GASV'
filenames=find_names(patname, getwd())
print(34)
print(filenames)
#开始统计
firstDF=T
for (bedname in filenames){
  print(bedname)
  #break}
  #bedname="SRR7707729.pindel.precision.overlap"
  group<-gsub('.overlap', '', bedname)
  #
  tryCatch({
	  bed=read.delim(bedname,header = F)
  }, warning = function(w){
    print("warning")
  }, error = function(e){
    print('the file is empty')
    print(bedname)
  },finally = {
    head(bed)
    #(head(bed))
    #
    colnames(bed)=c('chr1','start1','end1','cnvlen1','cnvtype1',
                    'chr2','start2','end2','cnvlen2','cnvtype2','cnvname','overlaplen')
    bed$cnvlen1=bed$end1-bed$start1
    head(bed)
    table(bed$cnvlen1)
    #统计每个bed中有多少CNV
    CNVnum=length(unique(paste(bed$chr1,bed$start1)))
    #bed$CNVtype=0
    ##real CNV name
    #unique(bed$CNVtype)
    head(bed)
    SUM1<-data.frame(matrix(0,2,5))
    #修改列名
    colnames(SUM1)<-c('group','cnvtype','positiveNUM','simuCNVNUM','positiveRATE')
    #对应关系
    #对应关系
    ###=================================================###
    ##CNVER的CNVTYPE
    bed$CNVtype1=0
    #del
    NAMEs=c('DEL')
    for (i in 1:length(NAMEs)){
      NAME=NAMEs[i]
      #赋值为del
      bed$CNVtype1[grep(NAME,bed$cnvtype1)]='del'
    }
    #dup
    NAMEs=c("DUP","INS")
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
    head(bed)
    ###=================================================###
    #专属部分——CNV名称结束
    head(bed)
    cnvtypes=c('del','dup')
    #生成相同的名称了
    for (i in 1:2){
      #i=1
      df1=bed[bed$CNVtype1==cnvtypes[i],]
      head(df1)
      #分母是检出该类型的所有CNV
      detectCNVNUM=length(unique(paste(df1$chr1,df1$start1,df1$end1)))
      #df1中全是该CNV类型
      b=df1[df1$CNVtype2==cnvtypes[i],]
      b
      #去重
      positiveNUM=length(unique(paste(b$chr1,b$start1,b$end1)))
      #df1中全是该CNV类型
      #匹配的CNV
      #找出CNVtype匹配的类型，为simuCNVtype
      #未去重的
      #去重之后的
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
  })
}
print(130)
rateSUM
rateSUM$group=gsub(".overlap|.sorted|.bam","",rateSUM$group)
rateSUM<-rateSUM[order(rateSUM$group),]
print(args[2])
#outname包括路径
outname=args[2]
write.csv(rateSUM, outname, row.names = F)
