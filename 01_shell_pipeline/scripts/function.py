import sqlite3
import multiprocessing
import xlrd
import xlwt
import datetime
from functools import partial
import openpyxl
import math

#-----base_function---
def rev_seq(seq): 
    #aim:将read2的序列反转成跟read1同一链上的方向与碱基排序
    #example:"ATCGTAA"-->"TTACGAT"
    comp_encoding = {'A':'T', 'T':'A','G':'C','C':'G','N':'N'}
    str_r2c = '' # reverse complement read2 sequence
    for n in range(0,len(seq)):
        str_r2c = str_r2c + comp_encoding[seq[n]] 
    str_r2rc = str_r2c[::-1] #逆序排列
    seq = str_r2rc # reverse read2 quality score
    return(seq)

def readr1(r1):
    lines = []
    lines += [r1.readline().strip()] + [r1.readline().strip()] + [r1.readline().strip()] + [r1.readline().strip()]
    if len(lines[0]) > 2:
        return lines
    else:
        return 'end of file'

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

def block_enrich(block,err,win,fq):
    #aim:判断固定序列是否存在于输入的fq之中
    #参数说明:block(固定序列),fq(以四单元list读如的fq),err(允许错配碱基数目),win(固定序列在fq序列上的位置范围)
    seq=fq[1]
    qua=fq[3]
    win_seq=seq[win[0]:win[1]]
    ts=block_in(block,win_seq,err)
    if ts!=-1:
        re_seq=seq
        re_qua=qua
        new_fq=[fq[0],re_seq,"+",re_qua]
    elif ts==-1:
        new_fq="none"
    return(new_fq)

def openr1(input_dir,filename):
    status = ""
    total_count = 0
    print(total_count)
    with open(f'{input_dir}/{filename}.fq','r') as r1:
            while True:
                if status == "complete":
                    break
                fq_block = []
                for n in range(0,1000000000):
                    line_block = readr1(r1)
                    print(line_block)
                    if line_block == "end of file":
                        print(f'end of {filename} file')
                        status = "complete"
                        break
                    elif len(line_block[1])>0 :
                        fq_block += [line_block]
                        total_count += 1
    return(fq_block)

def openr1_r2(input_dir,filename):
    status = ""
    total_count = 0
    with open(f'{input_dir}/{filename}_R1.fq', 'r') as r1: # read1 fastq file !!!后缀
            with open(f'{input_dir}/{filename}_R2.fq', 'r') as r2: # read2 fastq file !!!后缀
                while True:
                    if status == 'complete':
                        break
                    list_of_blocks = []
                    for n in range(0,100000000): #该循环用于将r1和r2每四行一一对应合并为一个字符串集,并计算所有的reads数目
                        lines_4x2 = read_paired_end(r1,r2)
                        if len(lines_4x2[1])>0 and len(lines_4x2[5])>0:
                            #lines_4x2[1]和lines_4x2[5]分别代表r1和r2的序列 
                            list_of_blocks += [lines_4x2]
                            total_count += 1
                        elif lines_4x2 == 'end of file':    
                            print('end of file')
                            status = 'complete'
                            break
    return(list_of_blocks)

def read_paired_end(r1,r2):
    lines_4x2 = []
    lines_4x2 += [r1.readline().strip()] + [r1.readline().strip()] + [r1.readline().strip()] + [r1.readline().strip()]
    lines_4x2 += [r2.readline().strip()] + [r2.readline().strip()] + [r2.readline().strip()] + [r2.readline().strip()]
    #r1.readline().strip()中readline()表示读取r1的下一行,strip()表示去除删除这行开头和结尾的空格和换行符，因此，代码的作用是读取文件对象 r1 的一行，并删除这行开头和结尾的空格和换行符，然后返回这一行内容的字符串
    #两行lines_4x2 +=...的作用是分别读取r1和r2的四行,组成字符串集合
    if len(lines_4x2[0]) > 2:
        return lines_4x2
    else:
        return 'end of file'

def write2r1(res_dir,filename,fq):
    with open(f'{res_dir}/{filename}.fq','a') as save_fastq:
        for bloc in fq:
            if bloc != "none":
                save_fastq.write(bloc[0] + '\n')
                save_fastq.write(bloc[1] + '\n')
                save_fastq.write(bloc[2] + '\n')
                save_fastq.write(bloc[3] + '\n')

def block_record(win,err,block_lis,fq):
    #function:根据输入的fq(四形式)以及block_lis,判断seq中是否存在block,输出格式为["fq[0]",1/0(1存在,0不存在),1/0...]
    #need_function:block_in
    seq=fq[1][int(win[0]):int(win[1])]
    block_in_seq=[]
    for i in range(0,len(block_lis)):
        block=block_lis[i]
        if block_in(block,seq,err)!=-1:
            block_in_seq.append(1)
        else:
            block_in_seq.append(0)
    block_in_seq.insert(0,fq[0])
    return(block_in_seq)

def block_record2xls(block_in_seqlis,block_lis,res_dir,xls_name):
    #function:将block_record的poolmap结果整理成xls表格
    block_lis.insert(0,"header")
    workbook = openpyxl.Workbook()
    sheet = workbook.active
    sheet.append(block_lis)
    for row in block_in_seqlis:
        sheet.append(row)
    workbook.save(f"{res_dir}/{xls_name}.xls")
    
    
#----function----------
def fastq_rev(lines_4X2): 
    #aim:将R2的fastq的list的序列反转成跟R1同一链的方向和剪辑排序
    #example:["@22","ACGTGTA","+","FIII:FF"] -> ['@22', 'TACACGT', '+', 'FF:IIIF']
    comp_encoding = {'A':'T', 'T':'A','G':'C','C':'G','N':'N'}
    str_r2c = '' # reverse complement read2 sequence
    for n in range(0,len(lines_4X2[1])):
        str_r2c = str_r2c + comp_encoding[lines_4X2[1][n]] 
    str_r2rc = str_r2c[::-1] #逆序排列
    lines_4X2[1] = str_r2rc
    lines_4X2[3] = lines_4X2[3][::-1]
    return(lines_4X2)

def merge_R1_frame(cut_score,frame,r1):
    #aim:将frame和r1进行合并
    #explanation:score(至少多少个碱基相同)
    #method:frame的第一个碱基为起点，r1的末尾碱基为起点进行比较,之后frame向前移动一个碱基,第一碱基与r1的倒二比较,第二与末尾比较……反复循环
    seq = r1[1]#获取r1的序列
    for n in range(cut_score-1,len(frame)):
        score=0
        s=0
        while s<=n:
            if frame[s]==seq[len(seq)-1-n+s]:
                score+=1
            s+=1
        frame_over=frame[0:s]
        r1_over=seq[len(seq)-1-n:]
        if score>cut_score-1 and score/len(frame_over)>0.75:#至少cut_score个碱基相同,且overlap相似率>75%,成功匹配
            new_seq = seq[0:len(seq)-1-n] + frame
            new_qua = r1[3][0:len(seq)-1-n] + "I"*len(frame) #用I补全质量行
            new_fastq = [r1[0],new_seq,"+",new_qua]
            break
    if n>=len(frame)-1:
        new_fastq="none"
    return(new_fastq)

def merge_frame_r2(cut_score,frame,r2):
    #aim:将frame与r2(反转后)进行合并,frame在前,r2在后
    #explanation:score(至少多少个碱基相同)
    #method:将frame与r2第1个碱基比较...，之后移动r2第1个碱基与frame第2个碱基比较...,循环往复
    seq=r2[1]
    pan=0
    for n in range(0,len(frame)):
        s=0
        score=0
        while s < len(frame)-n:
            if seq[s] == frame[s+n]:
                score+=1
            s+=1
        r2_only=seq[len(frame)-n:]
        frame_only=frame[0:n+1]
        r2_fra_over = seq[0:len(frame)-n]
        if score>cut_score-1 and score/len(r2_fra_over)>0.75:
            pan=1
            new_seq=frame+r2_only
            new_qua=(len(frame_only)-1)*"I"+r2[3]
            new_fastq=[r2[0],new_seq,"+",new_qua]
            break
    if pan==0:
        new_fastq="none"
    return(new_fastq)

def fastq_merge_r1_r2(r1r2,err):
    #aim:merge R1 and R2(fq) to merge_R
    #example
    workbook_qs = xlrd.open_workbook("quality_score.xls")
    worksheet_qs = workbook_qs.sheet_by_index(0)
    qs_encoding = {}
    for nrows in range(1,42):
        try:
            qs_code = int(worksheet_qs.cell_value(nrows,0))
        except:
            qs_code = worksheet_qs.cell_value(nrows,0)
        qs_encoding[str(qs_code)] = worksheet_qs.cell_value(nrows,2)
    mer_r1 = r1r2[0:4]
    mer_r2 = r1r2[4:8]
    rev_r2 = fastq_rev(mer_r2)
    seq_r1 = mer_r1[1]
    seq_r2 = rev_r2[1]
    qua_r1 = mer_r1[3]
    qua_r2 = rev_r2[3]
    pan = 0
    for n in range(0,len(seq_r1)+1):
        score = 0
        s = 0
        while s <= n:
            if seq_r2[s] == seq_r1[len(seq_r1)-1-n+s]:
                score+=1
            s+=1
    #-----get the fragment-----
        r2_over = seq_r2[0:s]
        r1_over = seq_r1[len(seq_r1)-1-n:]
        qua2_over = qua_r2[0:s]
        qua1_over = qua_r1[len(seq_r1)-1-n:]
        r2_seq_only = seq_r2[s:]
        r1_seq_only = seq_r1[0:len(seq_r1)-1-n]
        r2_qua_only = qua_r2[s:]
        r1_qua_only = qua_r1[0:len(seq_r1)-1-n]
    #---------------------------
        if score >= len(r2_over)-err and score/len(r2_over) > 0.75:
            pan = 1
            #-----over fragment-----
            base_loc = 0
            new_over_seq = ""
            new_over_qua = ""
            for base in r2_over:
                if qs_encoding[qua1_over[base_loc]] > qs_encoding[qua2_over[base_loc]]:
                    new_over_qua = new_over_qua + qua1_over[base_loc]
                    new_over_seq = new_over_seq + r1_over[base_loc]
                else:
                    new_over_qua = new_over_qua + qua2_over[base_loc]
                    new_over_seq = new_over_seq + r2_over[base_loc]
                base_loc += 1
        
            new_seq = r1_seq_only + new_over_seq + r2_seq_only
            new_qua = r1_qua_only + new_over_qua + r2_qua_only
            break
    if pan == 0:
        new_seq = "none"
        new_qua = "none"
    new_all = [r1r2[0],new_seq,"+",new_qua]
    return(new_all)

def seq2amacid(seq):
    #aim:将碱基翻译成氨基酸
    #example:"AAATTCGGTAGATTT" --> "KFGRF"
    matrix_cdo = {}
    workbook_cdo = xlrd.open_workbook("codon_table.xls")
    worksheet_cdo = workbook_cdo.sheet_by_index(0)
    for ncols_cdo in range(1,5):
        for nrows_cdo in range(1,17):
            cx = worksheet_cdo.cell_value(0,ncols_cdo)
            cy = worksheet_cdo.cell_value(nrows_cdo,0)
            cz = worksheet_cdo.cell_value(nrows_cdo,5)
            aa_cdo = worksheet_cdo.cell_value(nrows_cdo,ncols_cdo)
            codon_cdo = str(cy) + str(cx) + str(cz)
            matrix_cdo[codon_cdo] = str(aa_cdo)
    t=0
    orf=0
    if len(seq)%3 != 0:
        re="Error:Unequal number of bases to three."
    else:
        for base in seq:
            if base not in "ATCG":
                t=1
                re="Error:sequence contains illegal nucleotides"
        if t==0:
            data_aa = matrix_cdo.get(seq[0:3])
            while orf + 3< len(seq):
                orf=orf+3
                data_aa = data_aa + matrix_cdo.get(seq[orf:orf+3])
            re=data_aa
    return(re)
                    
def cut_seq(pre,behind,block,err,fq):
    #aim:将fq的固定序列及其前后几个碱基进行切除
    #参数说明:pre(固定序列的前n个碱基);behind(固定序列的后n个碱基);block(固定序列);err(允许错配碱基数目);fq(四单元)
    seq=fq[1]
    qua=fq[3]
    ts=block_in(block,seq,err)
    if ts!=-1:
        bef=ts-pre
        aft=ts+len(block)+behind
        if bef<0:
            bef=0
        if aft>len(seq):
            aft=len(seq)
        be_seq=seq[0:bef]
        af_seq=seq[aft:]
        be_qua=qua[0:bef]
        af_qua=qua[aft:]
        new_seq=be_seq+af_seq
        new_qua=be_qua+af_qua
        new_all=[fq[0],new_seq,"+",new_qua]
    else:
        new_all=["none"]
    return(new_all)

def extracr_seq_pre(pre,fq):
    #aim:取fq的前面pre个碱基序列
    #参数说明:pre(前n个碱基);fq(四单元)
    seq=fq[1]
    qua=fq[3]
    new_seq=seq[0:pre]
    new_qua=qua[0:pre]
    new_all=[fq[0],new_seq,"+",new_qua]
    return(new_all)

def extracr_seq_af(af,fq):
    #aim:取fq的后面af个碱基序列
    #参数说明:pre(前n个碱基);fq(四单元)
    seq=fq[1]
    qua=fq[3]
    new_seq=seq[-af:]
    new_qua=qua[-af:]
    new_all=[fq[0],new_seq,"+",new_qua]
    return(new_all)


def cut_seq_nomit(block_r1,err1,block_r2,err2,fq):
    #aim : 将fq的固定序列R1前及固定序列R2后进行切除
    #para : block(固定序列);err(允许错配碱基数目);fq(四单元),r1/r2(前面的序列和后面的序列)
    seq=fq[1]
    qua=fq[3]
    ts_r1=block_in(block_r1,seq,err1)
    ts_r2=block_in(block_r2,seq,err2)
    if ts_r1 != -1 and ts_r2 != -1:
        bef=ts_r1+len(block_r1)
        aft=ts_r2-1
        new_seq=seq[bef:aft]
        new_qua=qua[bef:aft]
        new_all=[fq[0],new_seq,"+",new_qua]
    else:
        new_all=["none"]
    return(new_all)

def cut_seq_after_nomit(block,err,fq):
    #aim : 将fq的固定序列及其后面的序列进行切除
    #para : block(固定序列);err(允许错配碱基数目);fq(四单元)
    seq=fq[1]
    qua=fq[3]
    end=block_in(block,seq,err)
    if end != -1:
        new_seq=seq[bef:end-1]
        new_qua=qua[bef:end-1]
        new_all=[fq[0],new_seq,"+",new_qua]
    else:
        new_all=["none"]
    return(new_all)

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


#----pipline-----
def block_desper(input_dir,filename,res_dir,block_lis,err,win,xls_name):
    #aim:将fq根据固定序列进行拆分
    #参数说明:input_dir(fq文件夹);filename(fq的name);res_dir(输出文件夹);block_lis(固定序列的集合);err(允许错配数目);win(扫描范围)
    #need:openr1,block_enrich,block_record,block_record2xls
    pool = multiprocessing.Pool()
    all_fq=openr1(input_dir,filename)
    for block_in in block_lis:
        tmp=pool.map(partial(block_enrich,block_in,err,win),all_fq)
        tmp=[x for x in tmp if x != ["none"]]
        with open(f'{res_dir}/{filename}_{block_in}.fq','a') as save_fastq:
            for bloc in tmp:
                if bloc != "none":
                    save_fastq.write(bloc[0] + '\n')
                    save_fastq.write(bloc[1] + '\n')
                    save_fastq.write(bloc[2] + '\n')
                    save_fastq.write(bloc[3] + '\n')
    in_read = pool.map(partial(block_record,win,err,block_lis),all_fq)
    block_record2xls(in_read,block_lis,res_dir,xls_name)

def block_cut(indir,outdir,pre,behind,block,err,name):
    #aim:将fq固定序列及其前后几个碱基进行切除
    #参数说明:indir(fq文件夹);outdir(输出文件夹);pre(固定序列前几个碱基);behind(固定序列后几个碱基);block(固定序列);err(允许错配碱基数目);name(fq名字)
    pool = multiprocessing.Pool()
    print(f'{indir}/{name}')
    ori_fq=openr1(indir,name)
    print(ori_fq[1])
    cut_fq=pool.map(partial(cut_seq,pre,behind,block,err),ori_fq)
    cut_fq=[x for x in cut_fq if x != ["none"]]        
    print(cut_fq[1])
    write2r1(outdir,name,cut_fq)

def block_cut_nomit(indir,inname,outdir,outname,threshold,block1,err1,block2,err2):
    #aim:将fq的R1固定序列及前面的序列切除，同时将fq的R2的固定序列及其后面碱基进行切除
    #参数说明:indir(fq文件夹);inname(输入文件名字);outdir(输出文件夹);outname(输出文件名字)
    #block1(固定序列R1);err1(允许错配碱基数目R1);block2(固定序列R2);err2(允许错配碱基数目R2);
    #name(输入和输出fq的名字),threshold(线程)
    pool = multiprocessing.Pool(threshold)
    print(f'{indir}/{inname}')
    ori_fq=openr1(indir,inname)
    print(ori_fq[1])
    cut_fq=pool.map(partial(cut_seq_nomit,block1,err1,block2,err2),ori_fq)
    cut_fq=[x for x in cut_fq if x != ["none"]]        
    print(cut_fq[1])
    write2r1(outdir,outname,cut_fq)


def merge_r1_fr_r2(frame,cut_score_r1,cut_score_r2,fastq):
    #--------------------
    #aim:将r1与frame合并,然后frame与r2合并,最后r1-fr和fr-r2合并
    r1=[fastq[0],fastq[1],fastq[2],fastq[3]]
    r2=[fastq[4],fastq[5],fastq[6],fastq[7]]
    r1_re = merge_R1_frame(cut_score_r1,frame,r1) #function!
    rev_r2 = fastq_rev(r2) #function!
    r2_re = merge_frame_r2(cut_score_r2,frame,r2) #function!
    if r1_re != "none" and r2_re != "none":
        res_seq = r1_re[1][0:(len(r1_re[1])-len(frame))] + r2_re[1]
        res_qua = r1_re[3][0:(len(r1_re[3])-len(frame))] + r2_re[3]
        res_fq = [r1_re[0],res_seq,"+",res_qua]
    else:
        res_fq="none"
    return(res_fq)

#------running--------
#------setting------
indir="/public/home/chenkai/shenlm/analysis/sk/241027_nss/res/1_cutt/tmp/"
outdir="/public/home/chenkai/shenlm/analysis/sk/241027_nss/res/1_cutt/"
sample="deumi" #the form is $sample".fq"
block1="CTAGTACGGGG"
threshold=10

#---beginning----
print('beiginning')
reads=openr1(indir,sample)
print(reads[1])
pool=multiprocessing.Pool(threshold)
la_seq=pool.map(partial(cut_seq_before_nomit,block1,1),reads)
write2r1(outdir,f'{sample}_cut',la_seq)








