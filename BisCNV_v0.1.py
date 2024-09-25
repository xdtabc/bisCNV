# -*- coding: utf-8 -*-
"""
Created on Fri Dem 10 9:28:14 2021

@author: Dantong Xu 
"""
#conda install code:
#conda install biopython numpy
#python code below:
from __future__ import division    #introduce the precise division 
import os, sys, re, argparse, random, subprocess #, gzip  #Warn: pysam only exists in Linux environment
#import pysam
import pandas as pd
import numpy as np
from Bio import SeqIO 
#from subprocess import Popen, PIPE
#os.chdir(r'D:\PRACTICE\sim_data\my_simu_fa')

##========================================================================================================================
##========================================================================================================================
##EXTERNAL INPUT AGRUMENTS 
parser = argparse.ArgumentParser(description='From FASTA to CNV\n')
parser.add_argument('--output_dir', 
                    help = 'The directory to generate output file\n')
#SIMULATE FASTA ARGUMENT
parser.add_argument('--SIMU_FASTA', action = 'store_true', default = False,    #THE FIRST MAIN ARGUMENT
                    help = 'CHOOSE TO GENERATE SIMULATED FASTA\n')
parser.add_argument('--input_fa', 
                    help = 'The reference fasta , whose extension must be .fa(or .fasta)\n')
parser.add_argument('--input_bed', 
                    help = 'The vcf, whose context include INSERTION and DELETION\n')

#SIMULATE FASTQ ARGUMENT
parser.add_argument('--SIMU_WGS_FASTQ', action = 'store_true', default = False,    #THE SECOND MAIN ARGUMENT
                    help = 'CHOOSE TO GENERATE WGS SIMULATED FASTQ\n')
parser.add_argument('--SIMU_WGBS_FASTQ', action = 'store_true', default = False,    #THE THIRD MAIN ARGUMENT
                    help = 'CHOOSE TO GENERATE WGBS SIMULATED FASTQ\n')
#The above 2 arguments control the 2nd and 3rd part of the program
parser.add_argument('--reference', 
                    help = 'The reference genome whose extension is .fa or .fasta , which must be input externally\n')
parser.add_argument('--seq_length', type=int, 
                    help = 'The length of sequence, please make sure it is between 0 bp and 300bp\n')
parser.add_argument('--conversion_rate', type = float, 
                    help = 'Use a uniform bisulfite conversionrate\n')
parser.add_argument('--CG_conversion', type = float,
                    help = 'Specify a CG-conversion rate between 0 and 100, this and the below argument "CH_conversion" is both used to generate the context specified bisulfite transformation\n')
parser.add_argument('--CH_conversion', type = float,
                    help = 'Specify a CH-conversion rate between 0 and 100\n')
parser.add_argument('--fixed_length_adapter', type = int, 
                    help = 'The length of adapter contamination must be shorter than the sequence length itself (and greater than 0)')
parser.add_argument('--variable_length_adapter', type = int,
                    help = 'The variable-length (insert-size depenedent) adapter contaminations\n')
parser.add_argument('--error_rate', type = float, default = 0.01, 
                    help = 'The error rate(%%), please make sure it is between 0 and 60%%, the dafault is 0.1(%%)\n')
parser.add_argument('--number_of_sequences', type = int, default = 1000000, 
                    help = 'The number of sequences to be generated\n')
parser.add_argument('--quality', type = int, default = 40, 
                    help = 'The quality values, please make sure it is in the range of (Phred) 2-40\n')
parser.add_argument('--non_directional', action = 'store_true', default = False, 
                    help = 'This will make non-directional BS-SEQ libraries\n')
parser.add_argument('--paired_end', action = 'store_true', default = False, 
                    help = 'This will make pair-end BS-SEQ libraries\n')
parser.add_argument('--minfrag', type = int, default = 70, 
                    help = 'The minimum fragment length, please make sure it is longer than 0 bp\n')
parser.add_argument('--maxfrag', type = int, default = 400, 
                    help = 'The maximum fragment length, please make sure it is shorter than 2000 bp\n')
parser.add_argument('--outfile_prefix', '-pre', 
                    help = 'The prefix of output file\'s name\n')

#DETECT CNV ARGUMENT
parser.add_argument('--DETECT_CNV', action = 'store_true', default = False,    #THE FOURTH MAIN ARGUMENT
                    help = 'CHOOSE TO DETECT CNV\n')
parser.add_argument('--fastq1', 
                    help = 'the name of input fastq1 file')
parser.add_argument('--fastq2',
                    help = 'The name of input fastq2 file')
#The form of input reference genome must be with a directory£¬such as /home/work/fasta/genome.fa
parser.add_argument('--genome',
                    help = 'The fasta file of the reference genome, end by .fa or .fasta')
#The form of output directory must end with '/', such as /home/work/output/
parser.add_argument('--chromosome', nargs='+', 
                    help = 'The chromosome which to detect CNV, the format needs to be "1, 2, 3,..., X, Y"')
parser.add_argument('--cnvnator_fa',
                    help = 'The path of storing CNVnator fasta file, the name of fasta needs to be "chr1.fa, chr2.fa, ..."')

#here the variables are global variable 
args = parser.parse_args()
output_dir = args.output_dir  

#SIMULATING FASTA's VARIABLE       #THE FIRST MAIN ARGUMENT
SIMU_FASTA = args.SIMU_FASTA
input_fa = args.input_fa
input_bed = args.input_bed

#SIMULATING FASTQ's VARIABLE
SIMU_WGS_FASTQ = args.SIMU_WGS_FASTQ                        #THE SECOND MAIN ARGUMENT
SIMU_WGBS_FASTQ = args.SIMU_WGBS_FASTQ                      #THE THIRD MAIN ARGUMENT
reference = args.reference                                  #this is essential
seq_length = args.seq_length                               #this is essential       
conversion_rate = args.conversion_rate                      #if choose 'SIMU_WGBS_FASTQ', only one from the two part: conversion_rate or (CG_conversion_rate, CH_conversion_rate), can and must be chosen 
CG_conversion_rate = args.CG_conversion                     
CH_conversion_rate = args.CH_conversion
fixed_length_adapter = args.fixed_length_adapter            #only one from the two part: fixed_length_adapter or variable_length_adapter, can and must be chosen 
variable_length_adapter = args.variable_length_adapter
error_rate = args.error_rate                                #this is optional
number_of_sequences = args.number_of_sequences              #this is optional
quality = args.quality                                      #this is optional
non_directional = args.non_directional                      #this is optional
paired_end = args.paired_end                                #this is optional
minfrag = args.minfrag                                      #this is optional
maxfrag = args.maxfrag                                      #this is optional
outfile_prefix = args.outfile_prefix                                        #this is optional

#DETECTING CNVS' VARIABLE          #THE SECOND MAIN ARGUMENT
DETECT_CNV = args.DETECT_CNV                                #this is essential
fastq1_name = args.fastq1                                   #this is essential
fastq2_name = args.fastq2                                   #this is essential
genome_name = args.genome                                   #this is essential
chromosome = args.chromosome                                #this is essential
cnvnator_fa = args.cnvnator_fa                              #this is essential

#防止测穿reads
if minfrag < seq_length:
    minfrag = seq_length

### OUTPUT DIRECTORY
if output_dir == None:     #if --output_dir isn't given, parser.parse_args().output_dir is None
    output_dir = os.getcwd()
    print("The output_dir hasn't been set, using the current folder:\n%s to output" %output_dir)

try:
    output_dir = os.path.abspath(output_dir)    #  making the path absolute
    os.chdir(output_dir)
except:
    print("Failed to open the output directory: %s!\n" %output_dir)
    sys.exit(0)
else:
    print("Output will be written into the directory: %s\n" %output_dir)

##========================================================================================================================
##========================================================================================================================
#SIMULATED FASTA FUNCTION

def create_the_output_filenames(input_bed):
    """
    function: create the output filenames
    parameter: input_fa, which stores the reference fasta and vcf
    """
    prefix = ''.join(input_bed.split('.')[:-1])
    simu_bed_name = prefix + '.simu.bed'
    simu_fa_name = prefix + '.simu.fa'
    return simu_bed_name, simu_fa_name
           
#this function is used to read the fasta file
def read_genome(fa_name):
    chromosomes = {}
    for seq_record in SeqIO.parse(fa_name, "fasta"):
        seq = seq_record.seq
        id = seq_record.id
        chromosomes[id] = seq   
    return chromosomes

#23:19 12/7/2021 i have too low power to die now
def generate_the_simulated_fasta(input_fa, input_bed):
    """
    function: generate the simulated fasta using the input reference fasta and vcf, containing the INSERSION, DUPLICATION and DELETION
              with the function 'find_input_file'
    parameter: workdir, which stores the reference fasta and vcf    
    """
    #input_fa='test_1108.fa'
    fa = read_genome(input_fa)    #fa is a dictionary of all chromosomes
    #open write_file
    simu_bed_name, simu_fa_name = create_the_output_filenames(input_bed) 
    #delete the old file
    if os.path.exists(simu_bed_name):
        os.remove(simu_bed_name)
    simu_bed = open(simu_bed_name,'w')  #open file just for 1 time!
    #use df.to_csv写入
    #open vcf file using pandas
    try:
        bed = pd.read_table(input_bed, header=None)
        #there must be 4 columns in the bed file
        bed.columns = ["chrom", "start", "length", "svtype"]   #rename the colnames 
    except:
        print("failed to read the bed file")
        sys.exit(0)
    #vcf=pd.read_table('test_1109_2.vcf', header=None)
    #按染色体划分，一条一条计算
    #修改DUP
    chrom_names = list(set(list(bed['chrom'])))
    for chrom_name in chrom_names:      
        print(chrom_name)  
        # break
        #读取原始的整条染色体序列
        #原本DUP的染色体
        seq_chrom = fa[chrom_name]    #here get the whole chromosomes sequence
        #尝试字符串分开相加，不用list
        simu_len = 0  #初始值为0；累加：每条染色体归零重新计算
    #每次计算1条染色体
        #x,y染色体时是否会有问题?
        try:
            chrom_name = int(chrom_name)
        except:
            chrom_name = chrom_name    #chrom_name是常染色体之外的序列，不用改变
        chr_f = bed.loc[bed['chrom'] == chrom_name]   #chr_f是从bed截取的来，有四列
        chr_f = chr_f.reset_index(drop = True)   #重新设置行索引
        #先处理所有的DUP
        for i in range(chr_f.shape[0]):  
            # i = 1
            chrom = chr_f.iloc[i,0]
            start = chr_f.iloc[i,1]  #start
            SV_length = chr_f.iloc[i,2]  #SV_length是原来bed的记录
            svtype = chr_f.iloc[i,3]
            #DUP不需要计算simu_st
            # simu_st = start + simu_len   #simu_st不需要累加 
            #先处理所有的DUP, 不考虑重叠的情况
            if svtype == '<DUP>':
                alt = seq_chrom[start: start + SV_length]   #提取DUP重复序列
                #a是除了本染色体之外剩余的染色体
                a = list(fa.keys()) 
                a.remove(str(chrom_name))
                #在剩余染色体个数中抽取三个随机数
                randlist = list(np.random.randint(0, len(a), size = 3))
                for rand_i in range(3):
                    #在剩余染色体a集合中按随机数抽取一条染色体
                    chrom_key = a[randlist[rand_i]]
                    chrom_fa = fa[chrom_key]     #fa不变 
                    rand = random.randint(1, len(chrom_fa))
                    simu_line = '\t'.join(list(str(x) for x in [chrom_key, rand, rand + SV_length, SV_length, svtype])) + '\n'
                    #写入simu_bed
                    simu_bed.write(simu_line)
                    print('alt', len(alt), 'SV_length', SV_length)
                    ####====================================###
                    chrom_fa = chrom_fa[: rand] + alt + chrom_fa[rand: ]
                    #改完要覆盖fa
                    fa[str(chrom_key)] = chrom_fa
    #修改INS、DEL
    for chrom_name in chrom_names:     
        print(chrom_name)  
        # break
        seq_chrom = fa[chrom_name]    #here get the whole chromosomes sequence
        chrom_len = len(seq_chrom)
        #尝试字符串分开相加，不用list
        simu_len = 0  #初始值为0；累加：每条染色体归零重新计算
        #每次计算1条染色体
        #x,y染色体时是否会有问题?
        try:
            chrom_name = int(chrom_name)
        except:
            chrom_name = chrom_name    #chrom_name是常染色体之外的序列，不用改变
        chr_f = bed.loc[bed['chrom'] == chrom_name]   #chr_f是从bed截取的来，有四列
        chr_f = chr_f.reset_index(drop = True)   #重新设置行索引
        #先处理所有的DUP
        for i in range(chr_f.shape[0]):  
            # i = 2
            chrom = chr_f.iloc[i,0]
            start = chr_f.iloc[i,1]  #start
            SV_length = chr_f.iloc[i,2]  #SV_length是原来bed的记录
            if chrom_len - start > abs(SV_length):
                SVLEN = SV_length  #INS、DEL的SVLEN与原来bed记录不变
                svtype = chr_f.iloc[i,3]
                # simu_st = 0 
                simu_st = start #+ simu_len   #simu_st不需要累加 
                #<CNV>当作<INS>处理
                if svtype == '<INS>' or svtype == '<CNV>':  #INSERTION: need 2 message(position, length), the alt sequences are generated randomly
                    #the alt sequence is generated in random
                    rand_list = list(np.random.randint(1, 5, size = SV_length))     #generate "SV_length" random numbers from 1 to 4
                    alt = []    
                    bases = {1:"A", 2:"T", 3:"C", 4:"G"}   #extract the bases randomly
                    for i in rand_list:
                        alt.append(bases[i])
                    alt = ''.join(alt)
                    #alt is recorded to be used for verification
                    #建立一个以fa为基准的bed文件
                    #改变的关键是simu_st，simu_ed是与simu_st绑定的
                    #SVLEN是计算截取区域用的长度：
                    simu_ed = simu_st + SV_length 
                    # simu_line = '\t'.join(list(str(x) for x in [chrom, simu_st, simu_ed, SV_length, svtype, alt, 'Y']))  #simu_st+1是必要的
                    simu_line = '\t'.join(list(str(x) for x in [chrom, simu_st, simu_ed, SV_length, svtype]))  
                    print("simu_line", simu_line)
                    #在INS_DEL_pos中记录START, END
                    # print('alt', len(alt), 'SV_length', SV_length)
                    #最后才写入bed文件，每条染色体写入一次
                    simu_bed.write(simu_line+'\n')
                    ####====================================###
                    seq_chrom = seq_chrom[: simu_st] + alt + seq_chrom[simu_st:]
                    #测试、检验
                elif svtype == '<DEL>':  #DELETION: need 2 message(position, length)
                    if SV_length < 0:
                        pass
                    else:
                        print('svtype: <DEL> is wrong with SV_length: %s' %SV_length)
                        sys.exit(0)
                    simu_ed = '-'
                    alt = '-'
                    # simu_line = '\t'.join(list(str(x) for x in [chrom, simu_st, simu_ed, SV_length, svtype, '-']))
                    simu_line = '\t'.join(list(str(x) for x in [chrom, simu_st, simu_ed, SV_length, svtype, 'Y']))
                    print(simu_line)
                    simu_bed.write(simu_line+'\n')
                    # print('SV_length', SV_length)
                    ####====================================###
                    #注意SV_length必须是负值
                    # seq_chrom = seq_chrom[: simu_st + SV_length] + seq_chrom[simu_st: ]         #定义SVLEN为缺失长度
                    seq_chrom = seq_chrom[: simu_st] + seq_chrom[simu_st - SV_length:]  
                else:
                    pass
                    #DUP不处理
                chrom_len = len(seq_chrom)
                #SVLEN的区别：INS、DEL与SV_length相等，DUP更复杂，如上所示
                #simu_len计算
                simu_len += SVLEN #simu_len需要累加    
                # count += 1
                #改完要覆盖fa
                fa[str(chrom_name)] = seq_chrom
            else:
                print("i", i, 'no')
                continue
    #检验
    # 每一条染色体独立地写入同一个fasta文件
    chrom_f = open(simu_fa_name,'w')
    # check_f = open('write_check', 'w')
    for chrom_name in fa.keys():
        print(chrom_name)
        chrom_fa = fa[chrom_name]
        #
        # check_f.write(chrom_name + '\t' + str(len(chrom_fa)) + '\n')
        header = '>' + str(chrom_name)
        chrom_f.write(header + '\n')
        # chr_len = 0
        for i in range(int(len(chrom_fa)/60)):
            newline = str(chrom_fa[60*i:60*(i + 1)])
            # chr_len += len(newline)
            chrom_f.write(newline + '\n')
        newline = str(chrom_fa[60*(i + 1):])
        # chr_len += len(newline)
        #记录每次写入染色体长度
        # check_f.write('write' + '\t' + str(chr_len) + '\n')
        chrom_f.write(newline + '\n')
    # check_f.close()
    #fasta写完
    simu_bed.close()
    chrom_f.close()    #chrom_f is the output simulated file
    print('simulating the fasta file has been done')

#generate_the_simulated_fasta(input_fa, input_bed)

##========================================================================================================================
##========================================================================================================================
#SIMULATED FASTQ FUNCTION
counting = {
    'non_directional': 0,
    'uniform_converted': 0,
    'uniform_total': 0,
    'CH_converted': 0,
    'CH_total': 0,
    'CG_converted': 0,
    'CG_total': 0,
    'small': 0,
    'total_base_count': 0,
    'total_errors_produced': 0
    }

#用于判断变量是否被定义的函数，不能用于外部传入参数的判断，这个函数的使用需谨慎
#now it can't be used for the argument
def isset(v): 
   try : 
     type (eval(v)) 
   except : 
     return 0  
   else : 
     return 1  

#produce_the_normal_random_number: this function generates random numbers with a normal distribution 
#loc(均值)，暂定为200，应该设置为输入参数；scale(标准差)，暂定为60，
def produce_the_normal_random_number():
    mu = variable_length_adapter   #this is the mean value of the length of variable adapters
    n = number_of_sequences
    sigma = 60   #sigma是方差, 暂时默认为60
    #n = number_of_sequences  #随机数的个数
    a = np.random.normal(loc = mu , scale= sigma**0.5, size = (1, n))    #生成随机正态分布数
    normal_random_numbers = [int(x) for x in list(a[0])]
    return normal_random_numbers   #normal_random_numbers相当于@gaussian
#normal_random_numbers为生成的具有n个正态分布随机数的集合

#produce_quality_values: this function generates quality values of sequence
def produce_quality_values(error_rate, seq_length):   #输出error_quality
    #var  #注意后面var的赋值
    if error_rate == 0:  #如果错误率为0，所有碱基质量值都是一致的(用同一个quality score)
        print("Error rate is 0, so that quality values will be constant--the Phred score of %s\n" %quality)
        quals = [quality]*seq_length   #储存序列所有碱基质量值，其个数与序列碱基个数相同，序列碱基个数是number_of_sequences
        #这里quals里每一个quality score都要转换为ASCII码
        quals = [change_phred_score_into_quality_string(qual) for qual in quals]   #Phred33体系，chr()转换为ASCII码
        no_error_quality = ''.join(quals)    #生成no_error碱基质量值的字符串
        return no_error_quality
    else:   
        print("Using a error rate of %s*0.01 to produce quality values\n" %error_rate)
        print("Calculating the slope of the error curve\n")
        #计算出error钟形曲线的斜率
        var = calculate_slope_of_the_error_rate_curve(error_rate)   #自定义函数：确定控制错误曲线的斜率
        # print("Error rates per bp will be modelled according to the formula:\n")
        # print("default base quality - 0.034286*position[bp] + 0.0009263*(position[bp]**2)) - 0.00001*(position[bp]**3)*%s)\n\n" %var)
        quals = [] 
        for x in range(1, seq_length+1):
            term1 = 0.034286*x
            term2 = 0.0009263*x**2
            term3 = 0.00001*var*x**3    
            rot_quality = quality - term1 + term2 - term3
            if rot_quality < 2:
                rot_quality = 2
            quals.append(rot_quality) 
        quals = [change_phred_score_into_quality_string(qual) for qual in quals]  #quals每个元素都被替换
        error_quality = ''.join(quals)
        return error_quality            

def change_phred_score_into_quality_string(qual):
    qual = int(qual)
    qual = chr(qual+33) 
    return qual

def calculate_slope_of_the_error_rate_curve(error_rate):  #输入参数是error_rate, 返回范围的上限？
    user_specified_error_rate = error_rate/100     #the argument error_rate means error_rate%
    lower_limit = 1  
    upper_limit = 1000000 
    old_upper_limit = upper_limit
    count = 0
    while True:
        count += 1
        var = lower_limit        
        var = upper_limit
        error_rate_upper_limit = calculate_mean_error_rate(var)
        half_distance = float('{:.4f}'.format((upper_limit - lower_limit)/2))
        if user_specified_error_rate <= error_rate_upper_limit:
            if error_rate_upper_limit - user_specified_error_rate <= 0.00000001:
                var = upper_limit
                break
            else:
                old_upper_limit = upper_limit
                upper_limit = half_distance + lower_limit
        elif user_specified_error_rate > error_rate_upper_limit:
            if user_specified_error_rate - error_rate_upper_limit <= 0.00000001:
                var = upper_limit
                break
            else:
                upper_limit = old_upper_limit                
                half_distance = float('{:.4f}'.format((upper_limit - lower_limit)/2))
                lower_limit = half_distance + lower_limit
        else:
            print("what else can there be? %s %s\n" %(user_specified_error_rate, half_distance))
            sys.exit(0)
        if var == lower_limit:    #结束迭代（不然循环不会终止）
            break
    return var
                
def calculate_mean_error_rate(var):     ##获取错误率的平均值        
    errors = []
    for x in range(1, seq_length+1):
        term1 = 0.034286*x
        term2 = 0.0009263*x**2
        term3 = 0.00001*var*x**3
        rot_quality = quality - term1 + term2 - term3
        if rot_quality < 2:
            rot_quality = 2
        error_rate = float('{:.4f}'.format(change_phred_score_into_error_probability(rot_quality))) 
        errors.append(error_rate)
    mean_error_rate = 0
    for rate in errors:
        mean_error_rate += rate
    mean_error_rate /= len(errors)
    return mean_error_rate

def change_phred_score_into_error_probability(qual):
    error_rate = 10**(-qual/10)
    return error_rate

def produce_sequencing_errors(seq, qual):   #seq与qual应该是等长的, used to generate the sequences containing sequencing errors
    total_base_count = 0
    total_errors_produced = 0
    bases = list(seq)
    quals = list(qual)
    if not len(bases) >= len(quals):
        print("The sequence length (\",%s,\") was shorter than the length of the quality string (\",%s,\")\n" %(len(bases), len(quals)))
        sys.exit(0)  #退出脚本
    for index in range(len(quals)):    #quals最后一个元素的索引   
        total_base_count += 1     #自增运算符：加1
        phred_score = change_quality_string_into_phred_score(quals[index])   #自定义函数：将quality字符（ASCII码）转换为phred分数
        error_rate = change_phred_score_into_error_probability(phred_score)  #自定义函数：将phred分数转换为碱基错误率      
        Random = random.randint(1,10000)   #Random会是1到10000的随机整数
        Random /= 10000   #Random会是0.0001到1之间的以0.0001为单位的随机小数
        if Random < error_rate:
            Random = random.randint(1, 3)   
            if bases[index] == 'A' or bases[index] == 'a':
                Random += 0
            elif bases[index] == 'T' or bases[index] == 't':
                Random += 1
            elif bases[index] == 'C' or bases[index] == 'c':
                Random += 2
            elif bases[index] == 'G' or bases[index] == 'g':
                Random += 3
            else:
                print("Here is the invalid base: \"%s[%s]\"\n" %(bases, index))
            Random %= 4  #求整除后的余数并赋值给变量本身
            #这个$random的取值完美避开了$bases[$index]原本的碱基，使得三种情况下替换碱基都不会与原本碱基一样
            DNA_bases = ['A', 'T', 'C', 'G']
            bases[index] = DNA_bases[Random]   #DNA_bases在开头定义了：DNA_bases = ['A', 'T', 'C', 'G']
            total_errors_produced += 1
    substituted_seq = ''.join(bases)
    counting['total_base_count'] += total_base_count
    counting['total_errors_produced'] += total_errors_produced
    return substituted_seq

def change_quality_string_into_phred_score(String):  #phred分数转换为错误可能性
    qual = ord(String) - 33
    return qual

def produce_genomic_sequences(total_genome_length, chromosomes):
    seq = ''
    coords = ''
    plus = 0
    minus = 0
    valid = 0
    while True:
        Random = random.randint(1, total_genome_length)  #用全基因组的范围产生随机数
        chromosome_length = 0   
        #chromosomes用get_genome_from_fasta生成
        for Chr in chromosomes:  #chromosomes 保证按染色体顺序排好
            chromosome_length += len(chromosomes[Chr])
            if Random + seq_length < chromosome_length:     #seq_length为外部参数
                if not len(chromosomes[Chr]) > chromosome_length - Random + seq_length:
                    break
                seq = chromosomes[Chr][chromosome_length - Random: chromosome_length - Random + seq_length].upper()
                if 'N' in seq:
                    break
                if len(seq) != seq_length:
                    break
                strand = random.randint(0,1)  
                if strand == 0:
                    seq = reverse_complement(seq)  #反向互补链
                    minus += 1
                else:
                    plus += 1
                valid += 1
                if plus == 1:
                    coords = Chr + str(chromosome_length - Random + 1) + '-' + str(chromosome_length - Random + seq_length - 1)  # shorten sequence by 1 on 3' end
                else:
                    coords = Chr + str(chromosome_length - Random + 1 + 1) + '-' + str(chromosome_length - Random + seq_length)  # shorten sequence by 1 on 5' end
        if valid >= 1:     
            break
    if plus == 1:
        return seq, 'plus', coords
    else:
        return seq, 'minus', coords

#this function is used to generate reverse complement strand
def reverse_complement(s):
    basecomplement = {
        "A":"T", "T":"A", "G":"C", "C":"G", "N":"N", "a":"t", "t":"a", "g":"c", "c":"g", "n":"n"
          }
    letters = list(s)
    letters = [basecomplement[base] for base in letters]
    letters =  ''.join(letters)
    return letters[::-1]  #[::-1]是倒序的

def get_genome_from_fasta(reference):
    total_genome_length = 0
    chromosomes = {}
    for seq_record in SeqIO.parse(reference,"fasta"):
        seq = seq_record.seq
        id = seq_record.id
        chromosomes[id] = seq   
    for value in chromosomes.values():
        total_genome_length += len(value)
    return total_genome_length, chromosomes

def produce_genomic_sequences_paired_end(minfrag, maxfrag, chromosomes, total_genome_length):
    fragment_length = random.randint(minfrag, maxfrag - 1)     #从minfrag到maxfrag-1的随机整数, fragment_length is useful!
    plus = 0
    minus = 0
    valid = 0
    while True:
        rand_num = random.randint(1, total_genome_length) 
        chrom_len = 0 
        for chrom in chromosomes:
            # break
            # print('chromosome:',chrom)
            chrom_len += len(chromosomes[chrom])
            if rand_num + fragment_length < chrom_len:
                #产生随机数时不一定轮到chromosomes中的第1条染色体
                seq_start = rand_num - (chrom_len - len(chromosomes[chrom]))
                seq1 = chromosomes[chrom][seq_start: seq_start + seq_length].upper()    #截取chromosomes[chrom]的seq_start位开始的seq_length个字符
                # seq_start = seq_start - 1
                seq2 = chromosomes[chrom][(seq_start + fragment_length - seq_length): (seq_start + fragment_length)].upper() 
                if 'N' in seq1 or 'N' in seq2:
                    break
                elif len(seq1) != seq_length or len(seq2) != seq_length:
                    break
                strand = random.randint(0,1)    
                if strand == 0:   
                    # seq_start = seq_start + 1
                    temp = reverse_complement(seq1) 
                    seq1 = reverse_complement(seq2) 
                    seq2 = temp
                    minus += 1
                    coords = chrom + ":" + str(seq_start + 2) + '-' + str(seq_start + fragment_length)    
                elif strand == 1:             #the seq1 and seq2 are at the same position as before
                    plus += 1  
                    coords = chrom + ":" + str(seq_start + 1) + '-' + str(seq_start + fragment_length - 1)
                valid += 1
                break 
        if valid >= 1: 
            break
    if plus == 1:
        return seq1, seq2, 'plus', coords
    else:   
        return seq1, seq2, 'minus', coords
    
def produce_non_directional_sequences(seq1):   #这个是用于single_end的
    RC = 0
    Random = random.randint(1, 100)
    if Random <= 50:
        RC = 1;
        if paired_end:   #paired_end是外部输入参数，只有是或否两个值
            pass
        else:
            seq1 = reverse_complement(seq1)
            counting['non_directional'] += 1
    if paired_end:   #paired_end同上
        return RC
    else:
        return seq1

#this function is used to generate bisulfite sequences uniformly
def changes_sequences_in_bisulfite_uniformly(seq):
    total_C_count = 0
    converted_C_count = 0
    bases = list(seq)
    for base in bases:
        if base == 'C':
            total_C_count += 1
            if conversion_rate==100:
                #不需要取随机数
                converted_C_count += 1
                bases = ['T' if i == base else i for i in bases]   
            else:
                Random = random.randint(1, 10000)/100     #生成从1开始到10000的随机整数，再除以100
                if Random <= conversion_rate:
                    converted_C_count += 1
                    bases = ['T' if i == base else i for i in bases]   #将bases里的'C'替换为'T'
    counting['uniform_converted'] += converted_C_count
    counting['uniform_total'] += total_C_count
    bisulfite_seq = ''.join(bases)
    return bisulfite_seq
                
    
#this function is used to generate bisulfite sequences with specifical context
def changes_sequences_in_bisulfite_specifically(seq):
    converted_CG_count = 0
    total_CG_count = 0
    converted_CH_count = 0
    total_CH_count = 0
    bases = list(seq)
    for index in range(len(bases)):
        if bases[index] == 'C':
            Random = random.randint(1, 10000)/100 
            if index + 1 >= len(bases):   #scalar@bases是@bases长度，即其元素个数
                total_CH_count += 1
                if Random <= CH_conversion_rate:    #CH_conversion_rate外部输入参数
                    converted_CH_count += 1
                    bases[index] = 'T'
            else:
                if bases[index + 1] == 'G':
                    total_CG_count += 1
                    if Random <= CG_conversion_rate:
                        bases[index] = 'T'
                        converted_CG_count += 1
                else:
                    total_CH_count += 1
                    if Random <= CH_conversion_rate:
                        bases[index] = 'T'
                        converted_CH_count += 1
    counting['CG_converted'] += converted_CG_count
    counting['CG_total'] += total_CG_count
    counting['CH_converted'] += converted_CH_count
    counting['CH_total'] += total_CH_count
    bisulfite_seq = ''.join(bases)
    bisulfite_seq = bisulfite_seq[0: 0 + (len(bisulfite_seq) - 1)]    #截取从0开始，length($bisulfite_seq)-1个字符
    return bisulfite_seq


def produce_the_fixed_length_adapter(seq, adapter_seq, fixed_length_adapter):    #generate fixed length adapter
    ### seq与adapter_seq是字符串    
    seq = list(seq)
    substitute_seq = list(adapter_seq[:fixed_length_adapter])    
    #fixed_length_adpater外部输入参数：取adapter_seq的前fixed_length_adapter位字符（碱基）
    #不需要在函数中输入，在脚本中作为全局参数
    seq = ''.join(seq[:-fixed_length_adapter] + substitute_seq)
    return seq


def produce_the_variable_length_adapter(seq_length, seq1, seq2, adapter_seq, normal_random_numbers):
    fragment = normal_random_numbers[0]
    seq_length = seq_length - 1   # the sequence has already been shortened by 1 bp again.
    if fragment < 0:
        fragment = 0
    if fragment < seq_length:
        counting['small'] += 1
        sub_length = seq_length - fragment
        substitute_seq = adapter_seq[0: 0 + sub_length]
        if paired_end:  #注意这是布尔运算值
            seq1 = ''.join(seq1[:-sub_length] + substitute_seq)     #对$seq1,从-$sub_length位开始的$sub_length个字符以$substitute_seq替换
            seq2 = ''.join(seq2[:-sub_length] + substitute_seq)
            if len(seq1) < seq_length:
                print("The sequence is now only \",%s ,\" bp long! $seq1\n" %len(seq1))
            if len(seq2) < seq_length:
                print("The sequence is now only \",%s ,\" bp long! $seq2\n" %len(seq2))
        else:
            seq1 = ''.join(seq1[:-sub_length] + substitute_seq)
            if len(seq1) < seq_length:
                print("The sequence is now only \",%s ,\" bp long! $seq1\n" %len(seq1))
    if paired_end:
        return seq1, seq2
    else:
        return seq1

#main function
def run_sequence_generation(reference, outfile_prefix):
    total_genome_length, chromosomes = get_genome_from_fasta(reference)   #chromosomes: this dictionary stored the names and sequences of the chromosomes, when running this function it will be generated of itself
    plus_strand_total = 0
    minus_strand_total = 0
    phred33_quality = produce_quality_values(error_rate, seq_length)  
    reference_name = reference.split('/')[-1]
    #区分两种情况：
    #一是不需要在外部参数输入结果文件名称，此时outfile_prefix由reference生成；
    #二是需要在外部参数输入结果文件名称，此时需要外部输入outfile_prefix参数，作为区别并行运行不同的结果，最后要将所有结果合并在一起
    if outfile_prefix == None:
        outfile_prefix = '.'.join(reference_name.split('.')[:-1])  
    else:
        pass  
    #路径分割符
    if '/' in output_dir:
        path_split = '/'
    else:
        path_split = '\\'
    if SIMU_WGS_FASTQ:
        print("producing WGS sequences")
        FQ1_NAME = output_dir + path_split + outfile_prefix + '.wgs_simulated_R1.fastq'
        FQ2_NAME = output_dir + path_split + outfile_prefix + '.wgs_simulated_R2.fastq'
        FQ_NAME = output_dir + path_split + outfile_prefix + '.wgs_simulated.fastq'
        if paired_end:        
            FASTQ_R1 = open(FQ1_NAME, 'w')    
            FASTQ_R2 = open(FQ2_NAME, 'w')
        else:
            FASTQ = open(FQ_NAME, 'w')       
    elif SIMU_WGBS_FASTQ:
        print("producing WGBS sequences")     
        FQ1_NAME = output_dir + path_split + outfile_prefix + '.wgbs_simulated_R1.fastq'
        FQ2_NAME = output_dir + path_split + outfile_prefix + '.wgbs_simulated_R2.fastq'
        FQ_NAME = output_dir + path_split + outfile_prefix + '.wgbs_simulated.fastq'
        if conversion_rate == None and None in (CG_conversion_rate, CH_conversion_rate):
            print("Only one from the two groups of the arguments: '--conversion_rate' or ('--CG_conversion_rate', '--CH_conversion_rate'), can and must be set in external input")
            sys.exit(0)
        if paired_end:        
            FASTQ_R1 = open(FQ1_NAME, 'w')    
            FASTQ_R2 = open(FQ2_NAME, 'w')
        else:
            FASTQ = open(FQ_NAME, 'w')   
    else:
        print("at least choose one of the two arguments: --SIMU_WGS_FASTQ, --SIMU_WGBS_FASTQ")
        sys.exit(0)
    adapter_seq = 'CAAGCAGAAGACGGCATACGAGATCGGTCTCGGCATTCCTGCTGAACCGCTCTTCCGATCT'
    while True:
        adapter_seq += adapter_seq
        if len(adapter_seq) > seq_length:
            break
    adapter_seq = reverse_complement(adapter_seq)     
    if variable_length_adapter:                
        normal_random_numbers = produce_the_normal_random_number()     #这个函数功能是生成具有n个正态分布随机数的集合: normal_random_numbers
    count = 0
    while True:
        if count % 100000 == 0:
            if not count == 0:
                print("%s sequences have been produced so far\n" %count)
        if paired_end:
            [seq1, seq2, strand, genomic_coords] = produce_genomic_sequences_paired_end(minfrag, maxfrag, chromosomes, total_genome_length) 
            len(seq1)
            len(seq2)
            if strand == 'plus':
                plus_strand_total += 1
            else:
                minus_strand_total += 1
            if SIMU_WGS_FASTQ:
                pass    #generate the WGS sequences without bisulfite conversion
            elif SIMU_WGBS_FASTQ:
                if CG_conversion_rate != None and CH_conversion_rate != None:    #判断是否存在这两个参数
                    seq1 = changes_sequences_in_bisulfite_specifically(seq1)
                    seq2 = changes_sequences_in_bisulfite_specifically(seq2)
                else:    #If there isn't 如果不存在就选择统一的转换（即CG、CH无差别转换）
                    seq1 = changes_sequences_in_bisulfite_uniformly(seq1)
                    seq2 = changes_sequences_in_bisulfite_uniformly(seq2)
            else:
                print("failed to generate WGS or WGBS sequences!")
                sys.exit(0)
            if non_directional:
                complemented = produce_non_directional_sequences(seq1)
                if complemented == 1:
                    seq2 = reverse_complement(seq2)   
                    temp = seq1
                    seq1 = seq2
                    seq2 = temp
                else:
                    seq2 = reverse_complement(seq2)  
            else:
                seq2 = reverse_complement(seq2)
        else:
            [seq1, strand, genomic_coords] = produce_genomic_sequences(total_genome_length, chromosomes)
            if strand == 'plus':
                plus_strand_total += 1
            else:
                minus_strand_total += 1
            if non_directional:
                seq1 = produce_non_directional_sequences(seq1)
            if SIMU_WGS_FASTQ:
                print("producing WGS sequences")     
            elif SIMU_WGBS_FASTQ:
                if CG_conversion_rate != None and CH_conversion_rate != None:
                    seq1 = changes_sequences_in_bisulfite_specifically(seq1)
                else:
                    if conversion_rate == 0:
                        seq1 = seq1[: len(seq1) - 1]    
                    else:
                        seq1 = changes_sequences_in_bisulfite_uniformly(seq1)
                print("producing WGBS sequences")      
        if fixed_length_adapter:
            seq1 = produce_the_fixed_length_adapter(seq1, adapter_seq, fixed_length_adapter)
            if paired_end:
                seq2 = produce_the_fixed_length_adapter(seq2, adapter_seq, fixed_length_adapter)
        elif variable_length_adapter:
            if paired_end:
                [seq1, seq2] = produce_the_variable_length_adapter(seq_length, seq1, seq2, adapter_seq, normal_random_numbers)
            else:
                [seq1] = produce_the_variable_length_adapter(seq_length, seq1, '', adapter_seq, normal_random_numbers)    
        if error_rate != 0:
            seq1 = produce_sequencing_errors(seq1, phred33_quality)
            if paired_end:
                seq2 = produce_sequencing_errors(seq2, phred33_quality)
        if paired_end:     
            # FASTQ_R1.write("@%s_%s_R1\n" %(count, genomic_coords))
            FASTQ_R1.write("@%s_%s\n" %(count, genomic_coords))
            FASTQ_R1.write("%s\n" %seq1)
            FASTQ_R1.write("+\n")
            FASTQ_R1.write("%s\n" %phred33_quality)
            # FASTQ_R2.write("@%s_%s_R2\n" %(count, genomic_coords))
            FASTQ_R2.write("@%s_%s\n" %(count, genomic_coords))
            FASTQ_R2.write("%s\n" %seq2)
            FASTQ_R2.write("+\n")
            FASTQ_R2.write("%s\n" %phred33_quality)
        else:
            FASTQ.write("@%s_%s\n" %(count, genomic_coords))
            FASTQ.write("%s\n" %seq1)
            FASTQ.write("+\n")
            FASTQ.write("%s\n" %phred33_quality)
        count += 1   
        if count == number_of_sequences:
            print("count: %s" %count)
            break
    print("Report:\n")
    print("-"*150,"\n")
    print("%s simulated sequences have been successfully produced in total (+ strand: %s\t - strand: %s)\n" %(count, plus_strand_total, minus_strand_total))
    if non_directional:
        percentage_non_directional = format(counting['non_directional']*100/count, ".2f")  #确保counting['non_directional']为数值型
        print("Sequences reverse complemented (non-directional library): %s, %s*0.01)\n\n" \
              %(counting['non_directional'], percentage_non_directional))
    if SIMU_WGS_FASTQ:
        print("wgs fastq file has been generated\n")
    elif SIMU_WGBS_FASTQ:
        if CG_conversion_rate != None and CH_conversion_rate != None:
            percentage_CG_conversion_rate = format(counting['CG_converted'] * 100 / counting['CG_total'], ".2f")
            print("Cytosines converted in CG-context: %s%%\n" %percentage_CG_conversion_rate)
            percentage_CH_conversion_rate = format(counting['CH_converted'] * 100 / counting['CH_total'], ".2f")
            print("Cytosines converted in CH-context: %s%%\n" %percentage_CH_conversion_rate)
        else:
            #print("CG_conversion_rate:", CG_conversion_rate, "CH_conversion_rate", CH_conversion_rate)
            if conversion_rate == 0:
                print("No bisfulfite conversion performed\n")
            else:
                percentage_uniform_C_conversion = format(counting['uniform_converted']*100/counting['uniform_total'], ".2f")
                print("here!")
                print("Cytosines bisulfite converted in any context: %s%%\n" %percentage_uniform_C_conversion)
    if variable_length_adapter:
        percent_small = format(counting['small']*100/count, ".2f")
        print("%s elements were smaller than the sequence length %s (%s%%) and had some adapter sequence produced\n" %(counting['small'], seq_length, percent_small))
    if fixed_length_adapter:
        print("The last %s of all sequences were replaced with adapter sequence\n" %fixed_length_adapter)
    if not error_rate == 0:    
        percentage_errors_produced = format(counting['total_errors_produced'] / counting['total_base_count'] * 100, ".2f")
        print("Random sequencing errors produced in total: %s (of %s bp in total) (percentage: %s)\n" %(counting['total_errors_produced'], counting['total_base_count'], percentage_errors_produced))
    #return   #this function has no return value, so that it can have no "return"
    #FASTQ_R1.close()    #no file.close() in case of error
    #FASTQ_R2.close()
    #FASTQ.close()

##========================================================================================================================
##========================================================================================================================
#DETECT CNV FUNCTION
#DETECT_CNV = True

#Prepare the name of alignment file made by bismark.   
def generate_bismark_filename(fastq1_name):                
    prefix = re.split('[./]', fastq1_name)[-2]   #fastq1_name, output_dir is external input argument
    pe_bam = output_dir + '/' + prefix + '_bismark_bt2_pe.bam'           
    print("prefix", prefix)                            
    sorted_pe_bam = output_dir + '/' + 'SORTED_' + prefix + '_bismark_bt2_pe.bam'
    ambig_pe_bam = output_dir + '/' + prefix + '_bismark_bt2_pe.ambig.bam'
    sorted_ambig_pe_bam = output_dir + '/' + 'SORTED_' + prefix + '_bismark_bt2_pe.ambig.bam' 
    merged_bam = output_dir + '/' + 'MERGED_' + prefix + '_bismark_bt2.bam'    
    return pe_bam, sorted_pe_bam, ambig_pe_bam, sorted_ambig_pe_bam, merged_bam

#Use bismark to align the WGBS data.						
#remember when the script is tested successfully, delete the # of 'bismark_genome_preparation' code
def deal_with_bismark(pe_bam, sorted_pe_bam, ambig_pe_bam, sorted_ambig_pe_bam, merged_bam):
    genome_folder = '/'.join(genome_name.split("/")[:-1])   #genome_name is external input argument
    stat_index = subprocess.getstatusoutput('bismark_genome_preparation --verbose %s' %genome_folder)     #genome_folder is external input argument
    if stat_index[0] == 0:     
        print('"bismark_genome_preparation --verbose %s" is done;' %genome_folder)
    else:
        print('"bismark_genome_preparation --verbose %s" is error;' %genome_folder)
        sys.exit(0)
    stat_align = subprocess.getstatusoutput('bismark --genome %s -1 %s -2 %s -o %s --ambiguous --ambig_bam' %(genome_folder, fastq1_name, fastq2_name, output_dir))
    if stat_align[0] == 0:
        print('"bismark --genome %s -1 %s -2 %s -o %s --ambiguous --ambig_bam" is done' %(genome_folder, fastq1_name, fastq2_name, output_dir))
    else:
        print('"bismark --genome %s -1 %s -2 %s -o %s --ambiguous --ambig_bam" is error' %(genome_folder, fastq1_name, fastq2_name, output_dir))
        sys.exit(0)
    print("pe_bam, sorted_pe_bam, ambig_pe_bam, sorted_ambig_pe_bam, merged_bam\n", pe_bam, sorted_pe_bam, ambig_pe_bam, sorted_ambig_pe_bam, merged_bam)
    os.system('ls pe_bam, sorted_pe_bam, ambig_pe_bam, sorted_ambig_pe_bam, merged_bam') 
    if os.system('samtools sort %s -O bam -T SORT -o %s && samtools index %s' %(pe_bam, sorted_pe_bam, sorted_pe_bam)) != 0:
        print("'samtools sort %s -O bam -T SORT -o %s && samtools index %s' is error" %(pe_bam, sorted_pe_bam, sorted_pe_bam))
        sys.exit(0)
    else:
        print("'samtools sort %s -O bam -T SORT -o %s && samtools index %s' is done" %(pe_bam, sorted_pe_bam, sorted_pe_bam))
    if os.system('samtools sort %s -O bam -T SORT -o %s && samtools index %s' %(ambig_pe_bam, sorted_ambig_pe_bam, sorted_ambig_pe_bam)):
        print("'samtools sort %s -O bam -T SORT -o %s && samtools index %s' is error" %(ambig_pe_bam, sorted_ambig_pe_bam, sorted_ambig_pe_bam))
        sys.exit(0)
    else:
        print("'samtools sort %s -O bam -T SORT -o %s && samtools index %s' is done" %(ambig_pe_bam, sorted_ambig_pe_bam, sorted_ambig_pe_bam))
    if os.system('samtools merge -f %s %s %s' %(merged_bam, sorted_pe_bam, sorted_ambig_pe_bam)) != 0:
        print('"samtools merge -f %s %s %s" is error' %(merged_bam, sorted_pe_bam, sorted_ambig_pe_bam))
        sys.exit(0)
    else:
        print('"samtools merge -f %s %s %s" is done' %(merged_bam, sorted_pe_bam, sorted_ambig_pe_bam))

#The function is used to open the sam or bam file.
def samOrBam(infile_name):                      #this is used in lots of functions
    if os.path.splitext(infile_name)[1] == '.sam':
        infile = open(infile_name,'r')
        return(infile)
    elif os.path.splitext(infile_name)[1] == '.bam':
        comBam = "samtools view -h %s" % infile_name
        p = subprocess.Popen(comBam, stdout = subprocess.PIPE, shell = True,universal_newlines = True)
        infile = p.stdout
        return(infile)

#The function is used to replace the string.
def replaceString(string, num, replace):       #this is used in change_base_in_sam()
    string2 = ''
    for i in range(len(string)):
        if i == num:
            string2 += replace
        else:
            string2 += string[i]
    return string2

#The function is used to replace the base in the reads of bam file.
def base_change(ref, target):                  #this is used in change_base_in_sam()
    #ref, target are str 
    for i in range(0, len(ref)):
        if (ref[i] == 'C') & (target[i] == 'T'):
            target = replaceString(target, i, 'C')
        elif (ref[i] == 'G') & (target[i] == 'A'):
            target = replaceString(target, i, 'A')
    return target

#Change the base of bam file.
def change_base_in_sam(infile_name):                           #this is used in RUN_THE_MAIN_PROCESS()
    #Using pysam to open fasta file.
    fa = pysam.FastaFile(genome_name)
    input_file = samOrBam(infile_name)
    changed_sam_name = '.'.join(infile_name.split(".")[:-1]) + '_changed.sam'
    output_sam = open(changed_sam_name, 'w')
    for line in input_file:
        if line.startswith('@'):
            output_sam.write(line)    
            continue                  #skip the header (the line start with '@')
        else:
            line = re.split(r'\t', line)
            chr = str(line[2])
            if line[5] == '*':
                continue
            elif chr.isalpha() and chr != 'X' and chr != 'Y':
                continue
            else:			
                cigar = re.split(r'([A-Z])', line[5])  
                if '' in cigar:
                    cigar.remove('')    
            seq = ''        
            start = 0
            ref_st = int(line[3])
            ref = ''              
            for i in range(0, len(cigar)):
                if (cigar[i] == 'S') | (cigar[i] == 'I'):
                    start = start + int(cigar[i-1])              
                elif cigar[i] == 'D':
                    ref_st = ref_st + int(cigar[i-1])            
                elif cigar[i] == 'M':  
                    ref = ref + fa.fetch(region=str(chr) + ":" + str(ref_st) + "-" + str(ref_st + int(cigar[i - 1]) - 1))  #找对应的参考基因组序列
                    ref_st = ref_st + int(cigar[i-1])   
                    seq = seq + line[9][start:(start + int(cigar[i-1]))]
                    start = start + int(cigar[i-1])     
            newseq = base_change(ref, seq)
            finseq = ''
            lstart = 0  
            nstart = 0  
            for i in range(0, len(cigar)):
                if (cigar[i] == 'S') | (cigar[i] == 'I'):  
                    finseq = finseq + line[9][lstart:(lstart + int(cigar[i-1]))]
                    lstart = lstart + int(cigar[i-1])            
                elif cigar[i] == 'M':    
                    finseq = finseq + newseq[nstart:(nstart + int(cigar[i-1]))] 
                    lstart = lstart + int(cigar[i-1])
                    nstart = nstart + int(cigar[i-1])
            line[9] = finseq     
            line = '\t'.join(line)          
            output_sam.write(line)     
    output_sam.close()
  
#12/13/2021 17:03
#using samtools to deal with sam file in the output of bam file
def handle_the_changed_samfile(merged_bam):                                                         #this is used in deal_with_bwa()
    changed_sam_name = '.'.join(merged_bam.split(".")[:-1]) + '_changed.sam'    #set the changed samfile name
    changed_bam_name = '.'.join(merged_bam.split(".")[:-1]) + '_changed.bam'
    changed_sorted_bam_name = '.'.join(merged_bam.split(".")[:-1]) + '_changed_sorted.bam'
    if os.system("samtools view -S %s -b -o %s" %(changed_sam_name, changed_bam_name)) != 0:
        print('The process of sam changed to bam makes wrong!')
        sys.exit(0)
    else:
        print('"The process of sam changed to bam" is done')
    if os.system("samtools sort -O bam -T SORT -n %s -o %s" %(changed_bam_name, changed_sorted_bam_name)) != 0:
        print("The process of bam sorted makes wrong!")
        sys.exit(0)
    else:
        print("The process of bam sorted is done")
         
def extract_reads(changed_fastq1_name, changed_fastq2_name, changed_sorted_bam_name):       #this is used in deal_with_bwa()
    changed_fastq1 = open(changed_fastq1_name, 'w')     #输出changed_fastq1文件
    changed_fastq2 = open(changed_fastq2_name, 'w')     #输出changed_fastq2文件
    changed_sorted_bam = samOrBam(changed_sorted_bam_name)          #输入bam文件为重新排序后的bam
    for l in changed_sorted_bam:
        if l[0] == '@':
            continue   #满足if条件，则跳过for循环剩余的语句，进入下一个循环
        else:
            line = re.split(r'\t', l)
            Flag = int(line[1])
            line1 = ('@' + line[0]).replace('_', ' ')      #bismark不需要添加name    
            line2 = line[9]
            line3 = '+'
            line4 = line[10]
            if (Flag == 99 or Flag == 83):
                fastq = [line1, line2, line3, line4]
                for row in fastq:
                    changed_fastq1.write(row + '\n')
            elif (Flag == 147 or Flag == 163):
                fastq = [line1, line2, line3, line4]
                for row in fastq:
                    changed_fastq2.write(row + '\n')
            else:   #the other Flag values is not accepted
                continue
    changed_fastq1.close()
    changed_fastq2.close()
    changed_sorted_bam.close()

def deal_with_bwa(output_dir, fastq1_name, fastq2_name, merged_bam):                         #this is used in RUN_THE_MAIN_PROCESS()
    #prepare bwa files'name and the changed_sorted_bam file
    changed_fastq1_name = output_dir + '/' + re.split("\.|\/", fastq1_name)[-2] + '_changed.fq'     
    changed_fastq2_name = output_dir + '/' + re.split("\.|\/", fastq2_name)[-2] + '_changed.fq'
    changed_sorted_bam_name = '.'.join(merged_bam.split(".")[:-1]) + '_changed_sorted.bam'   #generate the changed and sorted bamfile
    #Extract reads from changed bam to make fastq file.
    extract_reads(changed_fastq1_name, changed_fastq2_name, changed_sorted_bam_name)
    #bwa_process
    stat_bwa_index = subprocess.getstatusoutput('bwa index %s' %genome_name)
    if stat_bwa_index[0] == 0:
        print('"bwa index %s" is done' %genome_name)
    else:
        print('"bwa index %s" is error' %genome_name)
        exit(0)
    realign_sam = '.'.join(changed_sorted_bam_name.split(".")[:-1]) + '_realign.sam'
    realign_bam = '.'.join(changed_sorted_bam_name.split(".")[:-1]) + '_realign.bam'
    if os.system('bwa mem %s %s %s > %s && samtools view -b -S %s -o %s' %(genome_name, changed_fastq1_name, changed_fastq2_name, realign_sam, realign_sam, realign_bam)) != 0:
        print("'bwa mem %s %s %s > %s && samtools view -b -S %s -o %s' is error %(genome_name, changed_fastq1_name, changed_fastq2_name, realign_sam, realign_sam, realign_bam)")   
        sys.exit(0)
    else:
        print("'bwa mem %s %s %s > %s && samtools view -b -S %s -o %s' is done" %(genome_name, changed_fastq1_name, changed_fastq2_name, realign_sam, realign_sam, realign_bam))   
        
#using CNVnator to detect CNV with bam file
def deal_with_CNVnator(fastq1_name, output_dir, merged_bam):                   #this is used in RUN_THE_MAIN_PROCESS()
    changed_sorted_bam_name = '.'.join(merged_bam.split(".")[:-1]) + '_changed_sorted.bam'
    #prepare the required variables
    prefix = re.split('[./]', fastq1_name)[-2]
    realign_bam = '.'.join(changed_sorted_bam_name.split(".")[:-1]) + '_realign.bam'
    cnvnator_root = output_dir + '/' + prefix + '.root'
    cnvnator_cnv = output_dir + '/' + prefix + '.cnv'
    cnvnator_vcf = output_dir + '/' + prefix + '.vcf'
    chromosome_name = ' '.join([str(i) for i in chromosome])
    #run the cnvnator
    if os.system('ls %s' %realign_bam) != 0:
        print("here! wrong")
        sys.exit(0)
    print("""
            cnvnator_root, chromosome_name, realign_bam, cnvnator_fa, cnvnator_cnv
            """, cnvnator_root, chromosome_name, realign_bam, cnvnator_fa, cnvnator_cnv)
    if os.system('cnvnator -root %s -chrom %s -tree %s && \
                cnvnator -root %s -chrom %s -his 1000 -d %s && \
                cnvnator -root %s -stat 1000 && \
                cnvnator -root %s -partition 1000 && \
                cnvnator -root %s -chrom %s -call 1000 > %s' 
                %(cnvnator_root, chromosome_name, realign_bam, 
                  cnvnator_root, chromosome_name, cnvnator_fa, 
                  cnvnator_root, cnvnator_root, cnvnator_root, chromosome_name, cnvnator_cnv)) != 0:
        print("using cnvnator to detect CNV is error")
        sys.exit(0)
    os.system("export LC_ALL=C && chmod a+x %s" %(cnvnator_cnv))
    stat_cnvnator2VCF = subprocess.getstatusoutput('cnvnator2VCF.pl %s > %s' %(cnvnator_cnv, cnvnator_vcf))
    if stat_cnvnator2VCF[0] == 0:
        print('CNV calling has been done')
        return True
    else:
        print('cnvnator2VCF.py makes an error')
        sys.exit(0) 

##========================================================================================================================
##========================================================================================================================
#RUN_THE_MAIN_PROCESS()   #don't use function here, because there are some global variable need to be defined
if __name__ == '__main__':
    if SIMU_FASTA:
        #judge whether both '--input_fa' and '--input_bed' exist or not
        if None in (input_fa, input_bed):
            print("The arguments: '--input_fa', '--input_bed' must be set in external input")
            sys.exit(0)
        else:
            ##START TO SIMULATED FASTA
            generate_the_simulated_fasta(input_fa, input_bed)   #input_fa, input_bed两个参数输入
    elif SIMU_WGS_FASTQ or SIMU_WGBS_FASTQ:
        if seq_length == None:                     #其余的必需参数：例如二选一类型，太复杂，暂时不考虑
            print("The arguments: '--length' should be set in external input")
            sys.exit(0)
#        elif fixed_length_adapter == None and variable_length_adapter == None:
#            print("Only one from the two groups of the arguments: '--fixed_length_adapter' or '--variable_length_adapter', can and must be set in external input")
#            sys.exit(0)  #don't have to introduce adapter contaminants
        elif reference == None:
            print("The arguments: '--reference' must be set in external input")
            sys.exit(0)
        else:
            #this part is used both in SIMU_WGS_FASTQ and in SIMU_WGBS_FASTQ
            #this setting variables need to be put in the main function
            #counting is the dictionary of initial values
            #run WGS or WGBS sequences generation
            run_sequence_generation(reference, outfile_prefix)
    elif DETECT_CNV:
        if None in (fastq1_name, fastq2_name, genome_name, chromosome, cnvnator_fa,):
            print("The arguments: '--fastq1_name', '--fastq2_name', '--genome_name', '--chromosome', '--cnvnator_fa' must be set in external input")
            sys.exit(0)
        else:
            (pe_bam, sorted_pe_bam, ambig_pe_bam, sorted_ambig_pe_bam, merged_bam) = generate_bismark_filename(fastq1_name)
            deal_with_bismark(pe_bam, sorted_pe_bam, ambig_pe_bam, sorted_ambig_pe_bam, merged_bam)
            change_base_in_sam(merged_bam)   
            handle_the_changed_samfile(merged_bam)
            deal_with_bwa(output_dir, fastq1_name, fastq2_name, merged_bam)
            status = deal_with_CNVnator(fastq1_name, output_dir, merged_bam)
            if status:
                print("remove the temporary file")
                os.system('rm *.bam *.sam *.fq *.bai *.gz *.txt *.cnv')
    else:
        print('Only one of the four main arguments: "--SIMU_FASTA", "--SIMU_WGS_FASTQ", "--SIMU_WGBS_FASTQ", "--DETECT_CNV", must be chosen to input')
        sys.exit(0)
