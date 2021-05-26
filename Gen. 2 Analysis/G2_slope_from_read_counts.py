#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: MattDeJong
"""
################################
################################
#Read count files must be in same directory. Files named as "Gen_2_replicate[number]_[inducer concentration]mM.fasta"
#The code will return an array with mutant p-value, slope, y-intercept, slope standard deviation, and y-intercept standard deviation.
#DNA sequence based analysis
################################
################################

import numpy as np
import math
from scipy.stats import t as tstat

def translate(seq): 
       
    table = { 
        'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M', 
        'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T', 
        'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K', 
        'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',                  
        'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L', 
        'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P', 
        'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q', 
        'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R', 
        'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V', 
        'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A', 
        'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E', 
        'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G', 
        'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S', 
        'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L', 
        'TAC':'Y', 'TAT':'Y', 'TAA':'*', 'TAG':'*', 
        'TGC':'C', 'TGT':'C', 'TGA':'*', 'TGG':'W', 
    } 
    protein ="" 
    if len(seq)%3 == 0: 
        for i in range(0, len(seq), 3): 
            codon = seq[i:i + 3] 
            protein+= table[codon] 
    return protein

def IUPAC(seq):
    table = {'A':['A'],
             'C':['C'],
             'T':['T'],
             'G':['G'],
             'R':['A','G'],
             'Y':['C','T'],
             'S':['G','C'],
             'W':['A','T'],
             'K':['G','T'],
             'M':['A','C'],
             'B':['C','G','T'],
             'D':['A','G','T'],
             'H':['A','C','T'],
             'V':['A','C','G'],
             'N':['A','T','C','G'],
             ' ':['']}
    seq_lst = ['']
    for nt in seq:
        add_lst = table[nt]
        temp_lst = []
        for nt_add in add_lst:
            for seq_add in seq_lst:
                temp_lst += [seq_add + nt_add]
        seq_lst = temp_lst
    return(seq_lst)
       
def capitalize(seq):
    cap = {'a':'A',
        't':'T',
        'c':'C',
        'g':'G',
        'A':'A',
        'T':'T',
        'C':'C',
        'G':'G'}
    new_seq = ''
    for nt in seq:
        new_seq += cap[nt]
    return(new_seq)
        
WT_DNA = 'ATGGACAAACCACCGTACCTTCCGCGTCCTCGTCCGCCAAGACGTATTTACAACCGTTAATGA'
WT_AA = translate(WT_DNA)

def read_text(text_filename):
    text_file = open(text_filename, "r")
    seq = ''
    dic = {}
    count = 0
    for line in text_file:
        if line[0]=='>':
            name1 = 'A' + capitalize(seq[0:57]) +'AATGA'
            if name1 not in dic.keys():
                dic[name1] = count
            else:
                dic[name1] += count      
            seq_ID = line.rstrip()
            count = ''
            for i in range(len(seq_ID)-2,0,-1):
                if seq_ID[i] == '=':
                    break
                else:
                    count = seq_ID[i] + count
            count = int(count)
            seq = ''  
        else:
            seq += line.rstrip()
    text_file.close()
    return (dic)

sum_count_dic = {} #dictionary that contains read counts for each observed sequence 

for i,rep in enumerate([1,2,3]):
    for j,IPTG in enumerate(['0','0.1','0.2','0.5']):
        file = 'Gen_2_replicate'+str(rep)+'_'+str(IPTG)+'mM.fasta'
        dic = read_text(file)
        for key in dic.keys():
                count = dic[key]
                if key not in sum_count_dic.keys():
                    sum_count_dic[key] = np.ones([3,4])/10
                    sum_count_dic[key][i,j] = count
                else:
                    sum_count_dic[key][i,j] += count
             
filtered_sum_count_dic = {} #filtered read count dictionary ensuring mutants were observed at all uninduced conditions
for key in sum_count_dic.keys():
    counts = sum_count_dic[key]
    if counts[0,0] >= 1 and counts[1,0] >= 1 and counts[2,0] >= 1: #greater than or equal to one because sum count dictionary is intialized with read count = 1/10
        if np.sum(counts,axis=0)[0]>40: #Ensures mutants were observed at least 40 times total
            filtered_sum_count_dic[key] = sum_count_dic[key]

def total(dic):
    tot = np.zeros([3,4])
    for key in dic.keys():
        tot = tot + dic[key]
    return(tot)
    
total_obs = total(filtered_sum_count_dic) #provides an array with the total counts for each repcliate and condition

def frequency_calc(dic,total):
    dic_new={}
    for key in dic.keys():
            dic_new[key]=np.divide(dic[key],total)
    return(dic_new)

freq_dic = frequency_calc(filtered_sum_count_dic,total_obs) #frequency dictionary normalizing read counts per replicate and condition
     
def normalize(frequency_dic):
    new_dic = {}
    for key in frequency_dic.keys():
       arr = frequency_dic[key]
       new_dic[key] = np.zeros([3,4])
       for i in range(3):
           for j in range(4):
               new_dic[key][i,j] = arr[i,j]/arr[i,0]
    return(new_dic)

norm_dic = normalize(freq_dic) #Normalizes frequency by 0 mM condition for each replicate

Neg_lst = []   #negative control DNA sequence list for all negative controls assayed
for key in freq_dic.keys():
    if '*' in translate(key[0:57]):
        Neg_lst += [key]

def g_p(frequency_dic,Neg_lst):
    gp0= 9.85
    gp =np.zeros([3,4])
    for key in Neg_lst:
        for i in range(3):
            for j in range(4):
                gp[i,j] += gp0-math.log(frequency_dic[key][i,j],2)   
    gp = gp/len(Neg_lst)  
    return (gp)

gp_arr = g_p(freq_dic,Neg_lst) #bulk number of doubling times for each condition

def growth(normalized_dic,gp_array): 
    new_dic={}
    t=395
    for key in normalized_dic.keys():
        new_dic[key]=np.zeros([3,4])
        for i in range(3):
            for j in range(4):
                new_dic[key][i,j]=1/(t*math.log(math.e,2))*(math.log(normalized_dic[key][i,j],2)+gp_array[i,j])
    return(new_dic)
    
growth_rt_dic=growth(norm_dic,gp_arr) #growth rate for each variant at each induction condition and replicate

def potencymetric(dic):
    new_dic={}
    def slope(growth_rate_lst):
            b = growth_rate_lst[0]
            sum1 = 0
            sum3 = 0
            sum2 = 0
            x_lst = [0,0.1,0.2,0.5]
            for i,yi in enumerate(growth_rate_lst):
                xi = x_lst[i]
                sum1 += xi*yi
                sum2 += xi
                sum3 += xi**2
            return((2*sum1-2*b*sum2)/(2*sum3))
    for key in dic.keys():
        new_dic[key]=np.zeros(3)
        new_dic[key][0] = slope(dic[key][0,:])
        new_dic[key][1] = slope(dic[key][1,:])
        new_dic[key][2] = slope(dic[key][2,:])
    return(new_dic)

regression_dic=potencymetric(growth_rt_dic) #slope dictionary. Y-intercept was held constant.

def average_metric(dic):
    ave_dic={}
    for key in dic.keys():
        ave_dic[key]=(dic[key][0]+dic[key][1]+dic[key][2])/3   
    return(ave_dic)

average_dic=average_metric(regression_dic) #dictionary of average slope across three repcliates
    
def std_metric(regression_dic,average_dic): 
    std_dic={}
    for key in regression_dic.keys():
        slopes = regression_dic[key]
        slope_mean = average_dic[key]
        std_slope = 0
        for i in range(3):
            std_slope += (slopes[i]-slope_mean)**2
        std_slope = (std_slope/2)**0.5
        std_dic[key] = std_slope
    return(std_dic)

std_dic=std_metric(regression_dic,average_dic) #dictionary of slope standard deviation across three repcliates

def significance_test(av_dic,st_dic,WT_DNA,growth_rt_dic):
    new_dic={}
    def F_test(std1,std_WT): #returns boolean whether standard deviations are signifcantly different
        if std1>std_WT:
            if std1**2/std_WT**2>19:   
                return(True)
            else:
                return(False)
        else:
            if std_WT**2/std1**2>19:
                return(True)
            else:
                return(False)
    for key in av_dic.keys():
        slope = av_dic[key]
        WT_slope = av_dic[WT_DNA]
        std_dev = st_dic[key]
        WT_std_dev = st_dic[WT_DNA]
        
        if F_test(std_dev,WT_std_dev):  
            #t-test with unequal variance
            t=(slope-WT_slope)/(std_dev**2/3+WT_std_dev**2/3)**0.5
            dof=(std_dev**2/3+WT_std_dev**2/3)**2/((std_dev**2/3)**2/2+(WT_std_dev**2/3)**2/2)
            p_value=tstat.cdf(t,dof)  
            new_dic[key]=np.array([p_value,slope,growth_rt_dic[key][0,0],std_dev,0,1,t,dof]) 
        else:
            #t-test with equal variance
            new_dic[key]=np.zeros(3)
            s_pooled=((2*std_dev**2+2*WT_std_dev**2)/4)**0.5
            t=(slope-WT_slope)/s_pooled*(9/6)**0.5  
            dof=6-2
            p_value=tstat.cdf(t,dof)
            new_dic[key]=np.array([p_value,slope,growth_rt_dic[key][0,0],std_dev,0,0,t,dof])
    return(new_dic)

t_values=significance_test(average_dic,std_dic,WT_DNA,growth_rt_dic)    #1-tailed t-test relative to parental 

def p_value_array(dic):
    sort_lst=[]
    for i in dic.keys():
        sort_lst+=[(i,dic[i][0],dic[i][1],dic[i][2],dic[i][3],dic[i][4])]  #[name,p_value,aveslope,aveyint,stdevslope,stdevyint]
    return(sort_lst)

p_unsorted=p_value_array(t_values)  #list of p values

def sort(values):
    dtype=[('name', 'U63'), ('pvalue', 'float'),('average_slope', 'U30'),('average_yint', 'U30'),('std_dev_slope', 'U30'),('std_dev_yint', 'U30')]   #originally said ('std_dev_yint', float)
    a = np.array(values, dtype=dtype)      
    return(np.sort(a, order='pvalue'))

sorted_values=sort(p_unsorted)  #sort mutants by 1-tailed p-value

nameHandle = open('Gen_2_results_for_DNA.csv', 'w')  #writing .csv file
string = 'Sequence,p-value (1-tailed),Average slope (mM^-1min^-1),Average y-intercept (min^-1),Slope standard deviation (mM^-1min^-1),Y-intercept standard deviation (min^-1)'
nameHandle.write(string + '\n')
for i in range(len(sorted_values)):
    row = sorted_values[i]
    name = row[0]
    pval = str(row[1])
    slope = row[2]
    yint = row[3]
    slope_stdev = row[4]
    yint_stdev = row[5]
    string = name + ',' + pval  + ',' + slope  + ',' + yint  + ',' + slope_stdev  + ',' + yint_stdev
    nameHandle.write(string + '\n')
nameHandle.close()
    