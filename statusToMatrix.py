import numpy as np
import pandas as pd
import sys

def count2b1(num):
    c = 0
    while (c<num):
        num &= (num-1)
        c+=1
    return c

def statusToMatrix(names, output):
    content = {}
    #print(names,output)
    index = []
    for i, name in enumerate(names):
        index.append(name[:name.find('.')])
        with open(name) as f:
            lines = f.readlines()
        for line in lines:
            c = line.strip().split()
            c[1] = int('1'+c[1],2)
            if c[0] in content:
                content[c[0]][i] = c[1]
            else:
                content[c[0]] = {i:c[1]}
    
    l = len(names)
    s = np.zeros((l,l))
    c = np.zeros((l,l))
    for loc in content:
        value = content[loc]
        for i in value:
            for j in value:
                if i>j:
                    dis = count2b1(value[i]^value[j])
                    s[i,j]+=dis
                    s[j,i]+=dis
                    c[i,j]+=1
                    c[j,i]+=1
    s = s/c
    result = pd.DataFrame(data=s,index=index,columns=index)
    result.to_csv(output+'.csv',na_rep='NA')
if __name__=="__main__":
    #print(count2b1(7))
    output = sys.argv[1]
    names = sys.argv[2:]
    statusToMatrix(names,output)
