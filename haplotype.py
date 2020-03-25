import sys
import subprocess
import pysam
import os

def search(s,arr):
    l=0
    r=len(arr)-1
    #print(s,arr)
    while l<r:
        mid=(l+r)//2
        if arr[mid]<s:
            l = mid+1
            continue
        if arr[mid]>=s:
            r = mid-1
        #print(l,r)
    return (l+r)//2

def main():
    '''
    target_flanking_len=80
    rcTable={}
    rcTable{'A'}='T'
    rcTable{'T'}='A'
    rcTable{'G'}='C'
    rcTable{'C'}='G'
    rcTable{'N'}='N'
    rcTable{'R'}='Y'
    rcTable{'Y'}='R'
    rcTable{'M'}='K'
    rcTable{'K'}='M'
    rcTable{'S'}='S'
    rcTable{'W'}='W'
    '''
    strand_map={
            0:{'+':['CG','TG'],'-':['CG','CA']},
            1:{'-':['CG','TG'],'+':['CG','CA']}
            }
    #'+'(++,+-) or '-'(-+,--) means strand, +CG,-CG means methylated +C/-C; +TG/-CA means unmethylated +C/-C
    #position of C in CG in different strand do not change! Because bsmap transfer all reads into ++ and __ which all 
    #in the + strand. 


    cpgfile = sys.argv[1]
    bamfile = sys.argv[2]
    fn=bamfile
    cg_monitor = {}
    with open(cpgfile) as f:
        lines = f.readlines()
    dic={}
    os.system('rm '+fn+'.status')
    for line in lines:
        chr,s,e = line.strip().split()[:3]
        if chr in dic:
            dic[chr].append(int(s))
        else:
            dic[chr]=[int(s)]
        cg_monitor[chr+'_'+s] = [0,0]
    
    for chr in dic:
        dic[chr].sort()
    

    #finish build the binary look up table

    f=pysam.AlignmentFile(fn, "r")
    result=[]
    lastreads=[]
    dup_dic=set()
    for line in f:
        name = line.query_name
        flag = int((line.flag&(0x10))>0)
        chr = line.reference_name
        if chr=='chrM':
            continue
        s = line.reference_start
        reads = line.seq
        qua = line.query_qualities
        mismatch = int(line.get_tag('NM'))
        strand = line.get_tag('ZS') 
        #if strand[0]!='+' and strand[0]!='-' or strand[1]!='+' and strand[1]!='-':
        #    strand = temp[13][-2:]
        if not chr in dic: continue
        pos=search(int(s),dic[chr])
        e = s+len(reads)
        cpg_arr=[]
        for i in range(pos,pos+100):
            if i==len(dic[chr]) or  dic[chr][i]>e:
                break
            num=dic[chr][i]
            if num>=s and num<e:
                cpg_arr.append(num)
        #print(cpg_arr)
        #print(s)
        #print(reads)
        cpg_status=''
        containc=0
        containt=0
        cpg_pos=''
        #Here's a new method about strand specific haplotype
        strand_specific_table = strand_map[flag][strand[1]]
        seq_err=False
        for num in cpg_arr: 
            if num-s+1>=len(reads) or num-s-1<=1:
                continue
            cpg_in_reads=reads[num-s-1:num-s+1]
            mark = chr+'_'+str(num)
            #print(cg_monitor[mark])
            if cpg_in_reads==strand_specific_table[0]:
                if cg_monitor[mark][0]!=0 and cg_monitor[mark][0]!=cg_monitor[mark][1]:
                    seq_err=True
                    break
                cpg_status=cpg_status+'C'
                cpg_pos = cpg_pos+str(num)+','
                cg_monitor[mark][0]+=1 
                cg_monitor[mark][1]+=1
                containc=1
            else:
                if cpg_in_reads==strand_specific_table[1]:
                    if cg_monitor[mark][0]!=0 and cg_monitor[mark][0]!=cg_monitor[mark][1]:
                        seq_err=True
                        break
                    cpg_status=cpg_status+'T'
                    cpg_pos = cpg_pos+str(num)+','
                    cg_monitor[mark][0]+=0
                    cg_monitor[mark][1]+=1
                    containt=1
                    #Position should be different in +/- strand for cpg. I will deal with this later.
                else:
                    cpg_status=cpg_status+'N'
                    cpg_pos = cpg_pos+str(num)+','
                    #print(cpg_in_reads,strand_specific_table,num,s)
                    #seq_err=True
                    #break
        #print(cpg_status,seq_err)
        #print(temp,flag)
        #print(seq_err)
        if not cpg_status: continue
        if seq_err:
            num = str(num)
            while num in result[-1]:
                result.pop()
            continue
        result.append(str(chr)+'\t'+str(s)+'\t'+str(e-1)+'\t'+str(name)+'\t'+str(cpg_status)+'\t'+str(strand)+'\t'+str(cpg_pos[:-1])+'\t'+str(flag)+'\t'+str(mismatch)+'\n')
        '''
        if cpg_status:
            dup_s=chr+str(s)+str(e-1)+cpg_status+strand+cpg_pos[:-1]
            if dup_s in dup_dic:
                continue
            else:
                dup_dic.add(chr+str(s)+str(e-1)+cpg_status+strand+cpg_pos[:-1])
                lastreads.append(dup_s)
                if len(lastreads)>30:
                    remove_s = lastreads[0]
                    lastreads=lastreads[1:]
                    dup_dic.remove(remove_s)
            result.append(str(chr)+'\t'+str(s)+'\t'+str(e-1)+'\t'+str(name)+'\t'+str(cpg_status)+'\t'+str(strand)+'\t'+str(cpg_pos[:-1])+'\t'+str(flag)+'\t'+str(mismatch)+'\n')
        '''
        if len(result)>1000000:
            with open(fn+'.status','a') as ff:
                ff.writelines(result)
            result=[]
    f.close()
    with open(fn+'.status','a') as ff:
        ff.writelines(result)
        #print(fn,sum,mix,mix/sum)


if __name__=="__main__":
    #print(search(10542,[10469, 10471, 10484, 10489, 10493, 10497, 10525, 10542, 10563, 10571]))
    main()
