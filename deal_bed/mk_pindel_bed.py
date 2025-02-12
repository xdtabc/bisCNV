# -*- coding: utf-8 -*-
"""
Created on Sun Mar 27 23:35:34 2022

@author: 98127
"""

import os, sys
dirName=sys.argv[1]
os.chdir(dirName)


path1=os.getcwd()
def getFileName(path, ext):
    ''' 获取指定目录下的所有包含指定子串的文件名 '''
    names = []
    f_list = os.listdir(path)
    # print f_list
    for name in f_list:
        # 以.为分隔符
        name1=name.split('.')
        if ext in name1:
            names.append(name)
    return names
            
ext='vcf'  #vcf为指定子串
names=getFileName(path1, ext)
for name in names:
    pre=name.rstrip('.vcf')
    bedname=pre+'.pindel.bed'
    with open(name, 'r') as f, open(bedname, 'w') as out:  
        for line in f:
            #去除#开头的行（去除注释）
            if not line.startswith('#'):
                l=line.split('\t')
                chrom=l[0]
                start=l[1]
                a=l[7].split(';')
                #这个文件的注释是不工整的，因此用字符串查找（in）来赋值
                for b in a:
                    if 'END=' in b:
                        end=b.lstrip('END=')
                    elif 'SVLEN=' in b:
                        svlen=b.lstrip('SVLEN=')
                    elif 'SVTYPE=' in b:
                        svtype=b.lstrip('SVTYPE=')
                newl=[str(x) for x in [chrom, start, end, svlen, svtype]]
                newline='\t'.join(newl)+'\n'
                out.write(newline)
            else:
                continue
            