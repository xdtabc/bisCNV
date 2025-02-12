setwd('E:\\wgs_wgbs_file\\filter_cnvbed\\画图代码\\CNVlength\\input_data\\originfile')
dir()


patname='simuDATA'
#bedDIR=paste(baseDIR,'\\',patname,sep='')
#realDATA
outnames=c()
#当前目录下查找文件
for (name in dir()){
  #name="realDATA.breakdancer.recall.average.csv"
  a=strsplit(name, split = '[.]')[[1]]
  #break}
  if (patname %in% a ){
    #如果包含指定子字符串，返回T；否则返回F
    #添加到一个向量中
    outnames<-c(outnames, name)
  }else{
    next
  }
}
filenames=outnames
filenames


firstdf=T
for (filename in filenames){
  #filename=filenames[2]
  print(filename)
  df=read.csv(filename)
  #加软件名称
  if (firstdf){
    df_sum=df
    firstdf=F
  }else{
    df_sum=rbind(df_sum,df)
  }
}
head(df_sum)
df_sum$`100-1,000`=df_sum$X100.500+df_sum$X500.1000
df_sum$`1,000-10,000`=df_sum$X1000.5000+df_sum$X5000.1e4
df_sum=df_sum[,c(1:3,10:11,8:9)]
colnames(df_sum)[c(3,6:7)]=c('0-100','10,000-100,000','100,000-1,000,000')

#调整输入格式
newdf=data.frame(matrix(0,5,3))
colnames(newdf)=c('group','cnvlen','cnvnum')
newdf
head(newdf)
#df=df[,-2]
cnvlens=c('0-100','100-1,000','1,000-10,000',
          '10,000-100,000','100,000-1,000,000')
for (j in 1:length(cnvlens)){
  newdf[j,2]=cnvlens[j]
}
firstdf=T
for (i in 1:nrow(df)){
  newdf[,1]=df[i,1]
  for (j in 1:nrow(newdf)){
    newdf[j,3]=df[i,j+1]
  }
  if (firstdf){
    alldf=newdf
    firstdf=F
  }else{
    alldf=rbind(alldf,newdf)
  }
}
alldf
name='realDATA.average.CNVlength.sum.1.csv'
write.csv(alldf,name,row.names = F)

#以下为求平均值
#查找出相同比对软件和相同测试的四个样本
a=unique(df_sum$group)
a
#a=gsub('sherman_','',a)
a=gsub('sample|_1_|_2_|_3_|_4_','',a)
groupnames=unique(a)
#此时groupnames即为每个组合的名称
groupnames
#计算平均值
df=df_sum
colnames(df)[4]='CNVNUM'
df
firstDF<-T

for (j in 1:length(groupnames)){
  #j=1
  groupname=groupnames[j]
  df1=df[grep(groupname,df$group),]
  df1
  #i=1
  SUM1<-data.frame(matrix(0,1,9))
  colnames(SUM1)=c('group','CNVnum','0-100','100-500','500-1000','1000-5000',
                   '5000-1e4','1e4-1e5','1e5-1e6')
  SUM1[1,1]=groupname
  for (j in 2:ncol(df1)){
    SUM1[1,j]=mean(df1[,j])
  }
  if (firstDF){
    allSUM<-SUM1
    #关闭开关
    firstDF<-F
  }else{
    allSUM<-rbind(allSUM, SUM1)
  }
}
allSUM
allSUM<-allSUM[order(allSUM$group),]
setwd('E:\\wgs_wgbs_file\\filter_cnvbed\\average_csv')
name='simuDATA.average.CNVlength.csv'
write.csv(allSUM, name,row.names = F)
