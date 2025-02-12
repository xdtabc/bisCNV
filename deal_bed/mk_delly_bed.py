# -*- coding: utf-8 -*-
"""
Created on Mon Mar 28 10:14:53 2022

@author: 98127
"""

import os, sys
vcfDIR=sys.argv[1]
os.chdir(vcfDIR)
names=os.listdir()
names

filepath=os.getcwd()
suffix='.vcf'
def getFileName(filepath, suffix):
    ''' 
    获取指定目录下的所有指定后缀的文件名 
    filepath参数为指定目录
    suffix参数为指定后缀
    '''
    namesout = []
    f_list = os.listdir(filepath)
    # print f_list
    for i in f_list:
        # os.path.splitext():分离文件名与扩展名
        if os.path.splitext(i)[1] == suffix:
            namesout.append(i)
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


dectorNAME='delly'
names=getFileName(filepath, suffix)

for name in names:
    pre=name.split('.')[0]
    bedname=pre+'.'+dectorNAME+'.bed'
    with open(name, 'r') as f, open(bedname, 'w') as out:  
        for line in f:
            if line.startswith('#'):
                continue
            else:
                l=line.split('\t')
                chrom=l[0]
                start=l[1]
                a=l[7].split(';')
                #这个文件的注释是不工整的，因此用字符串查找（in）来赋值
                #更新：in不能保证有且只有匹配字符串
                end, svlen, svtype='', '', ''
                for b in a:
                    if is_number(b.lstrip('END=')):                    
                        end=b.lstrip('END=')
                    elif 'SVTYPE=' in b:
                        svtype=b.lstrip('SVTYPE=')
                svlen=int(end)-int(start)
                newl=[str(x) for x in [chrom, start, end, svlen, svtype]]
                newline='\t'.join(newl)+'\n'
                out.write(newline)