import numpy as np
import pandas as pd
import sys
import os

def count2b1(num):
    c = 0
    while (c<num):
        num &= (num-1)
        c+=1
    return c

def statusToMatrix(names, output):
    content = {}
    #print(names,output)
    if os.path.exists(output+'.region.txt'):
        os.remove(output+'.region.txt')
    with open(names) as f:
        lines = f.readlines()
    for i in range(len(lines)):
        lines[i] = lines[i].strip()
    names = lines
    index = []
    length = {}
    for i, name in enumerate(names):
        index.append(name[:name.find('.')])
        with open(name) as f:
            lines = f.readlines()
        for line in lines:
            c = line.strip().split()
            s = c[1]
            c[1] = int('1'+c[1],2)
            if c[0] in content:
                content[c[0]][i] = c[1]
                length[c[0]] = len(s)
            else:
                content[c[0]] = {i:c[1]}
                length[c[0]] = len(s)
    
    l = len(names)
    s = np.zeros((l,l))
    c = np.zeros((l,l))
    result=[]
    for loc in content:
        value = content[loc]
        for i in value:
            for j in value:
                if i>j:
                    dis = count2b1(value[i]^value[j])
                    s[i,j]+=dis/length[loc]
                    s[j,i]+=dis/length[loc]
                    c[i,j]+=1
                    c[j,i]+=1
                    cc = loc.strip().split('-')
                    temp = [cc[0],cc[1],str(i),str(j),bin(value[i])[2:],bin(value[j])[2:],str(dis/length[loc])]
                    result.append('\t'.join(temp)+'\n')
        if len(result)>1000000:
            with open(output+'.region.txt','a') as ff:
                ff.writelines(result)
    if len(result)>0:
        with open(output+'.region.txt','a') as ff:
            ff.writelines(result)
    s = s/c
    result = pd.DataFrame(data=s,index=index,columns=index)
    result.to_csv(output+'.csv',na_rep='NA')

if __name__=="__main__":
    #print(count2b1(7))
    output = sys.argv[1]
    names = sys.argv[2]
    statusToMatrix(names,output)
