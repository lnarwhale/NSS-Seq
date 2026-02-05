#!/usr/bin/python
import string
import argparse
import multiprocessing
from functools import partial
import math

parser=argparse.ArgumentParser(description="smart3_cut of the smart-seq3-data")
parser.add_argument('-i','--input',help='the location of fastq before cut')
parser.add_argument('-o','--output',help='the location of output')
parser.add_argument('-n','--name',help='the name of fastq')
args=parser.parse_args()

threshold=20
indir=args.input
outdir=args.output
sample=args.name

def readr1(r1):
    lines = []
    lines += [r1.readline().strip()] + [r1.readline().strip()] + [r1.readline().strip()] + [r1.readline().strip()]
    if len(lines[0]) > 2:
        return lines
    else:
        return 'end of file'

def openr1(input_dir,filename):
    status = ""
    total_count = 0
    #print(total_count)
    with open(f'{input_dir}/{filename}.fq','r') as r1:
            while True:
                if status == "complete":
                    break
                fq_block = []
                for n in range(0,1000000000):
                    line_block = readr1(r1)
                    #print(line_block)
                    if line_block == "end of file":
                        print(f'end of {filename} file')
                        status = "complete"
                        break
                    elif len(line_block[1])>0 :
                        fq_block += [line_block]
                        total_count += 1
    return(fq_block)

def write2r1(res_dir,filename,fq):
    with open(f'{res_dir}/{filename}.fq','a') as save_fastq:
        for bloc in fq:
            if bloc != "none":
                save_fastq.write(bloc[0] + '\n')
                save_fastq.write(bloc[1] + '\n')
                save_fastq.write(bloc[2] + '\n')
                save_fastq.write(bloc[3] + '\n')

def block_in(block,seq,err):
    #aim:判断block是否在错误允许(err)中存在于seq中，是的话返回匹配的第一个碱基的位置,否则返回-1
    #method:n是seq的碱基逐个递进,s是block上的碱基逐个递进,在n的循环中循环s,从而判断从第n个碱基开始的len(block)个碱基是否与block相同
    pan=0
    for n in range(0,len(seq)-len(block)+1):
        score=0
        s=0
        while s<len(block):
            if block[s]==seq[s+n]:
                score+=1
            s+=1
        if score>=len(block)-err:
            pan=1
            re=n
            break
    if pan==0:
        re=-1
    return(re)

def cut_seq_before_nomit(block,err,fq):
    #aim : 将fq的固定序列及其前面的序列进行切除
    #para : block(固定序列);err(允许错配碱基数目);fq(四单元)
    seq=fq[1]
    qua=fq[3]
    end1=block_in(block,seq,err)
    end=end1+len(block)
    if end1 != -1:
        new_seq=seq[end:]
        new_qua=qua[end:]
        new_all=[fq[0],new_seq,"+",new_qua]
    else:
        new_all=["none"]
    return(new_all)

'''find block1'''
fil=open("{0}/{1}.fq".format(indir,sample),"r")
time=1
n1=0
n2=0
block_1="CTAGTACGGGG"
block_2="TCGCCTTAGGG"
for pan in fil:
    loc1=pan.find(block_1)
    loc2=pan.find(block_2)
    if loc1 != -1:
        n1=n1+1
    if loc2 != -1:
        n2=n2+1
    time=time+1
    if time>80:
        break
if n1>n2:
    block1=block_1
else :
    block1=block_2
fil.close()

reads=openr1(indir,sample)
pool=multiprocessing.Pool(threshold)
la_seq=pool.map(partial(cut_seq_before_nomit,block1,1),reads)
write2r1(outdir,f'{sample}',la_seq)

