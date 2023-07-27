#!/usr/bin/env python3
#####miRNA Nucleotide bias analysis and length distribution analysis
import sys

def read_fa(fa):
    with open(fa) as f:
        for line in f.readlines():
            if not line.startswith('>'):
                line=line.strip()
                yield line

def prefer_nucl(inputfa,outputPath):
    prefer_nucl_dic={}
    for line in read_fa(inputfa):
            mi_length=len(line)
            for i in range(mi_length):
                position=i+1
                prefer_nucl_dic.setdefault(position,[]).append(line[i])
    with open(f'{outputPath}/prefer_nucl.csv','w') as w:
        for key in prefer_nucl_dic.keys():
            U_prefer_percent=prefer_nucl_dic[key].count('U')/len(prefer_nucl_dic[key])*100
            A_prefer_percent=prefer_nucl_dic[key].count('A')/len(prefer_nucl_dic[key])*100
            C_prefer_percent=prefer_nucl_dic[key].count('C')/len(prefer_nucl_dic[key])*100
            G_prefer_percent=prefer_nucl_dic[key].count('G')/len(prefer_nucl_dic[key])*100
            w.write(str(key)+','+str(U_prefer_percent)+','+'U'+'\n')
            w.write(str(key)+','+str(A_prefer_percent)+','+'A'+'\n')
            w.write(str(key)+','+str(C_prefer_percent)+','+'C'+'\n')
            w.write(str(key)+','+str(G_prefer_percent)+','+'G'+'\n')

def first_nucl_bias(inputfa,outputPath):
    first_nucldic={}
    for line in read_fa(inputfa):
        mi_length=len(line)
        first_nucldic.setdefault(mi_length,[]).append(line)
    with open(f'{outputPath}/first_nucl_bias.csv','w') as w:
        for key in first_nucldic.keys():
            first_content=[]
            for i in first_nucldic[key]:
                first_content.append(i[0])
            U_prefer_percent=first_content.count('U')/len(first_nucldic[key])*100
            A_prefer_percent=first_content.count('A')/len(first_nucldic[key])*100
            C_prefer_percent=first_content.count('C')/len(first_nucldic[key])*100
            G_prefer_percent=first_content.count('G')/len(first_nucldic[key])*100
            w.write(str(key)+','+str(U_prefer_percent)+','+'U'+'\n')
            w.write(str(key)+','+str(A_prefer_percent)+','+'A'+'\n')
            w.write(str(key)+','+str(C_prefer_percent)+','+'C'+'\n')
            w.write(str(key)+','+str(G_prefer_percent)+','+'G'+'\n')                                


def length_distri(inputfa,outputPath):
    length_list=[]
    for line in read_fa(inputfa):
        miRNA_length=len(line)
        length_list.append(miRNA_length)

    max_value=max(length_list)
    min_value=min(length_list)
    with open(f'{outputPath}/length_Distribution.csv','w') as w:
        for i in range(min_value,max_value+1):
            Length_distribution=length_list.count(i)
            w.write(str(i)+','+str(Length_distribution)+'\n')

if __name__=='__main__':
    args=sys.argv
    prefer_nucl(args[1],args[2])
    length_distri(args[1],args[2])
    first_nucl_bias(args[1],args[2])