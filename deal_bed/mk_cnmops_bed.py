# -*- coding: utf-8 -*-
"""
Created on Mon Mar 28 10:14:53 2022

@author: 98127
"""

import os, sys
dirName=sys.argv[1]
os.chdir(dirName)


filepath=os.getcwd()
suffix='.out'
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

names=getFileName(filepath, suffix)

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


dectorNAME='cn_mops'
for name in names:
    pre=name.split('.')[0]
    bedname=pre+'.'+dectorNAME+'.bed'
    with open(name, 'r') as f, open(bedname, 'w') as out:  
        for line in f:
            if line.startswith('seqnames'):
                continue
            else:
                l=line.split('\t')
                line='\t'.join(l[:4]+[l[-1]])
                out.write(line)
