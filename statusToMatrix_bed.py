import numpy as np
import pandas as pd
import sys
import argparse
from datetime import datetime


def detect_ham_dis(status_arr, target, pos_info, mats,  keeplarge=True):
    if keeplarge:
        f = target+1
    else:
        f = 0
    chr,cpg_pos,adjacent_cgs = pos_info
    mat, frac_number = mats
    pre_cg, next_cg = adjacent_cgs
    target_len = len(status_arr[target])
    result = [-1] * len(status_arr)
    outlines = []
    for i in range(f,len(status_arr)):
        if i!=target and status_arr[i]!='':
            min_length = min(len(status_arr[i]), target_len)
            result[i] = 0
            for j in range(min_length):
                if status_arr[target][j]!=status_arr[i][j]:
                    result[i]+=1
            result[i] /= float(min_length)
            # print(i,target,result)
            # print(mat)
            # print(frac_number)
            mat[i,target]+=result[i]
            mat[target,i]+=result[i]
            frac_number[i,target]+=1
            frac_number[target,i]+=1
            s, e = int(cpg_pos[-1*min_length:][0]), int(cpg_pos[-1*min_length:][-1])
            if next_cg!='None':
                distance_next = int(next_cg)-e
            else:
                distance_next = 'None'

            if min_length>=len(cpg_pos):
                if pre_cg!='None':
                    distance_prev = s - int(pre_cg)
                else:
                    distance_prev = 'None'
            else:
                distance_prev = s - int(cpg_pos[-1*min_length-1])
            output = [chr, ','.join(cpg_pos[-1*min_length:]),str(distance_prev),str(distance_next), str(target), str(i), status_arr[target][:min_length], status_arr[i][:min_length], str(result[i])]
            outlines.append('\t'.join(output)+'\n')
    outf.writelines(outlines)
    

def detect_all(status_arr, pos_info, mats):
    #chr,cpg_pos = pos_info
    for i in range(len(status_arr)):
        if status_arr[i]!='':
            detect_ham_dis(status_arr, i, pos_info, mats)



def statusToMatrix(cpg_file, names, desert_thres=1000):
    print('Start: ',str(datetime.now()))

    with open(cpg_file) as f:
        lines = f.readlines()
    
    chrs, cpgs = [], []
    temp = []
    for line in lines:
        if 'chrM' in line: continue
        c = line.strip().split()
        if len(chrs)==0 or c[0]!=chrs[-1]:
            if len(chrs)!=0:
                cpgs.append(temp)
            chrs.append(c[0])
            temp = []
        temp.append(str(int(c[1])-1))
    cpgs.append(temp)
   
    
    identified_cg = {}
    with open(names) as f:
        ns = f.readlines()
    for i,name in enumerate(ns):
        name = name.strip()
        with open(name) as f:
            lines = f.readlines()
        for line in lines[1:]:
            if 'chrM' in line: continue
            c = line.strip().split()
            chr,start,end,ratio = c[:4]
            if ratio!='0' and ratio!='1': continue
            start = start
            if chr in identified_cg:
                if start in identified_cg[chr]:
                    identified_cg[chr][start].append([i,ratio])
                else:
                    identified_cg[chr][start] = [[i,ratio]]
            else:
                identified_cg[chr] = {start:[[i,ratio]]}
    print('Finish data loading: ',str(datetime.now()))
    
    names = ns
    for i in range(len(names)):
        names[i] = names[i].strip()
    mat, frac_number = np.zeros((len(names),len(names))), np.zeros((len(names),len(names)))
    # print(identified_cg) 
    status = ['']*len(names)
    for i, chr in enumerate(chrs):
        cpg_pos = cpgs[i]
        for j, pos in enumerate(cpg_pos):
            if j-1>=0 and int(pos)-int(cpg_pos[j-1])>desert_thres:
                max_len = max(list(map(lambda x:len(x),status)))
                if max_len==0: continue
                if j-max_len-1<0:
                    prev_cg = 'None'
                else:
                    prev_cg = cpg_pos[j-max_len-1]
                detect_all(status, [chr,cpg_pos[j-max_len:j],(prev_cg,cpg_pos[j])], (mat, frac_number))
                status = ['']*len(names)
            update_set = set()
            if pos in identified_cg[chr]:
                for (target, ratio) in identified_cg[chr][pos]:
                    #print(target,ratio)
                    status[target]+=ratio
                    update_set.add(target)
                for k in range(len(names)):
                    if k not in update_set and status[k]!='':
                        if j-len(status[k])-1<0:
                            prev_cg = 'None'
                        else:
                            prev_cg = cpg_pos[j-len(status[k])-1]
                        detect_ham_dis(status, k, [chr,cpg_pos[j-len(status[k]):j],(prev_cg,cpg_pos[j])], (mat, frac_number), False)
                        status[k]=''
                #for k in range(len(names)):
                #    if k not in update_set:
                #        status[k]=''
            else:
                max_len = max(list(map(lambda x:len(x),status)))
                if max_len==0: continue
                if j-max_len-1<0:
                    prev_cg = 'None'
                else:
                    prev_cg = cpg_pos[j-max_len-1]
                detect_all(status, [chr,cpg_pos[j-max_len:j],(prev_cg,cpg_pos[j])],(mat, frac_number))
                status = ['']*len(names)

        max_len = max(list(map(lambda x:len(x),status)))
        if max_len==0: continue
        detect_all(status, (chr,cpg_pos[-1*max_len:],(cpg_pos[-1*max_len-1],'None')),(mat, frac_number))
        status = ['']*len(names)

    matrixFormation(mat,frac_number, names)
    print('Finish: ',str(datetime.now()))


def matrixFormation(mat, frac_number, names):
    result = pd.DataFrame(data=mat/frac_number,index=names,columns=names)
    result.to_csv(output+'.matrix.csv',na_rep='NA')

output=''

#mat = None
#frac_number = None

if __name__=="__main__":
    #print(count2b1(7))
    parser = argparse.ArgumentParser(description='Pairwise overlapped region for scWGBS data.')
    parser.add_argument('-cg','--cpgfile', help='cpgfile')
    parser.add_argument('-n','--names', help='scWGBS mcall G.bed file name list')
    parser.add_argument('-d','--desert_threshold', help='Maximum length between two adjacent CpG',type=int,default=1000)
    parser.add_argument('-o','--output',help="output file prefix")
    args = parser.parse_args() 
    output = args.output
    outf = open(output,'w')
    header = ['chr','cg_position','distance_to_prev_cg','distance_to_next_cg','order1','order2','status1','status2','HamDis']
    outf.writelines(['\t'.join(header)+'\n'])
    statusToMatrix(args.cpgfile, args.names, args.desert_threshold)
    outf.close()
