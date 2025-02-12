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
#outDIR='E:\\wgs_wgbs_file\\filter_cnvbed\\real_overlap\\precision'
#0324
#realDATA
#indir="/BIGDATA2/scau_xlyuan_1/xdt/WGS_WGBS_callingCNV/pipeline/CNV_result/CNVpre_rec/precision"
args <- commandArgs(trailingOnly = T)
print(args[1])
indir=args[1]
#indir="E:\\wgs_wgbs_file\\filter_cnvbed\\DRAW_CODE\\GRADUATION\\CNVpre_rec\\pig\\precision"
setwd(indir)
dir()
#setwd(outDIR)
patname='BreakDancer'
filenames=find_names(patname, getwd())
filenames
#开始统计
firstDF=T
for (bedname in filenames){
  print(bedname)
  tryCatch({
    bed=read.delim(bedname,header = F)
  }, warning = function(w){
    print("warning")
  }, error = function(e){
    print('the file is empty')
    print(bedname)
  },finally = {
    head(bed)
    #odd_rownum=seq(1,nrow(bed),2)
    #偶数行是重叠个数
    #even_rownum=seq(2,nrow(bed),2)
    #bed<-cbind(bed[odd_rownum,], bed[even_rownum,2])
    head(bed)
    #
    group<-gsub('.overlap', '', bedname)
    group<-gsub('precision', '', group)
    colnames(bed)=c('chr1','start1','end1','cnvlen1','cnvtype1',
                    'chr2','start2','end2','cnvlen2','cnvtype2','cnvname','overlaplen')
    head(bed)
    summary(bed$cnvlen1)
    group<-gsub('.precision.overlap', '', bedname)
    #
  #  colnames(bed)=c('chr1','start1','end1','cnvlen1','cnvtype1',
  #                  'chr2','start2','end2','cnvlen2','cnvtype2','overlaplen')
    #统计每个bed中有多少simuCNV
    simuCNVnum=length(unique(paste(bed$chr1,bed$start1)))
    a=sort(unique(bed$cnvtype1))
    print(a)   #}
    table(bed$cnvtype1)
    sort(unique(bed$cnvtype1))
    sort(unique(bed$cnvtype2))
    #重新命名
    #bed$CNVtype<-0
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
    cnvtypes=c('del','dup')
    for (i in 1:2){
      #i=1
      #detectCNV中所有符合的CNV类型
      df1=bed[bed$CNVtype1==cnvtypes[i],]
      df1
      #去重了，计算分母
      detectCNVNUM=length(unique(paste(df1$chr1,df1$start1,df1$end1)))
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
  })
  if (firstDF){
    rateSUM<-SUM1
    #关闭开关
    firstDF<-F
  }else{
    rateSUM<-rbind(rateSUM, SUM1)
  }
}
rateSUM
rateSUM$group=gsub(".overlap|.sorted|.bam","",rateSUM$group)
rateSUM<-rateSUM[order(rateSUM$group),]
print(args[2])
#outname包括路径
outname=args[2]
write.csv(rateSUM, outname, row.names = F)

