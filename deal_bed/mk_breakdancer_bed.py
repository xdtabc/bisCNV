# -*- coding: utf-8 -*-
"""
Created on Mon Mar 28 10:14:53 2022

@author: 98127
"""

import os, sys
dirName=sys.argv[1]
os.chdir(dirName)
os.listdir()

def getFileName(filepath, suffix):
    ''' 
    获取指定目录下的所有含有指定字符串的文件名 
    filepath参数为指定目录
    suffix参数为指定后缀
    '''
    namesout = []
    f_list = os.listdir(filepath)
    # print f_list
    for fname in f_list:
        # os.path.splitext():分离文件名与扩展名
        if suffix in fname:
            namesout.append(fname)
    return namesout

def is_number(s):
    '''
    Parameters
    ----------
    s : TYPE
        DESCRIPTION.

    Returns
    -------
    bool
        DESCRIPTION.

    '''
    try:
        float(s)
        return True
    except ValueError:
        pass
    try:
        import unicodedata
        unicodedata.numeric(s)
        return True
    except (TypeError, ValueError):
        pass
    return False

filepath=os.getcwd()
suffix='breakdancer'
names=getFileName(filepath, suffix)
dectorNAME='breakdancer'


for name in names:
    pre='.'.join(name.split('.')[:-1])
    bedname=pre+'.bed'
    with open(name, 'r') as f, open(bedname, 'w') as out:  
        for line in f:
            if line.startswith('#'):
                continue
            else:
            #两个位置都记下来
                l=line.split('\t')
                chrom=l[0]
                svtype=l[6]
                svlen=abs(int(l[7]))
                #有两个起点
                start1=l[1]
                end1=int(start1)+svlen
                start2=l[4]
                end2=int(start2)+svlen
                #两个位置都写下来
                newl1=[str(x) for x in [chrom, start1, end1, svlen, svtype]]
                newl2=[str(x) for x in [chrom, start2, end2, svlen, svtype]]
                newline1='\t'.join(newl1)+'\n'
                out.write(newline1)
                newline2='\t'.join(newl2)+'\n'
                out.write(newline2)
                
                
          
##检验是否存在end<start的情况——这是错误情况
#name='Sample1.bwa.breakdancer.bed'
#pre='.'.join(name.split('.')[:-1])
#bedname=pre+'.'+'new'+'.bed'
#with open(name, 'r') as f, open(bedname, 'w') as out:  
#    for line in f:                
#        l=line.split()
#        if int(l[1]) < int(l[2]):
#            out.write(line)
#        else:
#            print(line, 'wrong')
                

  
                
  