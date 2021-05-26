#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: MattDeJong
"""

################################
################################
#Read count files must be in same directory. Files named as "Gen_1_replicate[number]_[inducer concentration]mM.tsv"
#The code will return an array with mutant p-value, slope, y-intercept, slope standard deviation, and y-intercept standard deviation.
#DNA sequence based analysis
################################
################################

import csv
import numpy as np
import statistics as stat
import math 
from scipy.stats import t as tstat

def total():   
    total=15*[0]
    file_name_lst = []
    reps = [1,2,3]
    conditions = ['t0','0mM','0.15mM','0.30mM','0.50mM']
    for rep in reps:
        for condition in conditions:
            file_name_lst += ['Gen_1_replicate'+str(rep)+'_'+condition]
    for i in range(1,16):    #indices 1 through 15 for the files
        # set our filename
        input_tsv = file_name_lst[i-1] + '.tsv'
        with open(input_tsv) as fd:
            rd = csv.reader(fd, delimiter="\t", quotechar='"')
            for row in rd:
                if row[5]=='Counts':
                    pass
                else:
                    total[i-1]+=int(row[5])
    return(total)

total = total() #provides a list with the total counts for each repcliate and condition

def sum_count():
    dic={}
    file_name_lst = []
    reps = [1,2,3]
    conditions = ['t0','0mM','0.15mM','0.30mM','0.50mM']
    for rep in reps:
        for condition in conditions:
            file_name_lst += ['Gen_1_replicate'+str(rep)+'_'+condition]
    for i in range(1,16):
        input_tsv = file_name_lst[i-1] + '.tsv'
        with open(input_tsv) as fd:
            rd = csv.reader(fd, delimiter="\t", quotechar='"')
            for row in rd:
                if row[5]=='Counts':
                    pass
                else:
                    if int(row[5])>1:
                        if (row[3]) not in dic.keys():
                            dic[row[3]]=np.ones(15)  # set pseudo count equal to one
                            dic[row[3]][i-1]=int(row[5])  
                        else:
                            dic[row[3]][i-1]+=int(row[5])
    return(dic)

sum_count_dic = sum_count()  #dictionary that contains read counts for each observed sequence 

def quality_check(dic):
    dic_new={}
    for key in dic.keys():
        count=0
        if dic[key][1]>1:  #greater than one because sum count dictionary is intialized with read count = 1
            count+=1
        if dic[key][6]>1:
            count+=1
        if dic[key][11]>1:
            count+=1
        if count>2:
            dic_new[key]=dic[key]
    return(dic_new)

sum_count_dic = quality_check(sum_count_dic)  #filtered read count dictionary ensuring mutants were observed at all uninduced conditions

def frequency_calc(dic,total):
    dic_new={}
    for key in dic.keys():
        if np.sum(dic[key])>100:  #Ensures mutants were observed at least 100 times total
            dic_new[key]=np.zeros(15)
            for i in range(15):
                dic_new[key][i]=dic[key][i]/total[i]
    return(dic_new)

freq_dic = frequency_calc(sum_count_dic,total)  #frequency dictionary normalizing read counts per replicate and condition

def normalize_15(dic):  
    dic_new={}
    for key in dic.keys():
        dic_new[key]=np.zeros(15)
        for i in range(0,5):
            dic_new[key][i]=dic[key][i]/dic[key][1]
        for i in range(5,10):
            dic_new[key][i]=dic[key][i]/dic[key][6]
        for i in range(10,15):
            dic_new[key][i]=dic[key][i]/dic[key][11]
    return(dic_new)
      
normalized_dic = normalize_15(freq_dic)  #Normalizes frequency by 0 mM condition for each replicate
    
def g_p(dic):
    a=['ATGGATTGATAGCCGTATTTACCACGTCCCCGCCCCCCACGTCGCATTTACAACCGC','ATGGATTAACCATGATATTTACCACGTCCCCGCCCCCCACGTCGCATTTACAACCGC','ATGGATAAGCCATAATAATTACCACGTCCCCGCCCCCCACGTCGCATTTACAACCGC','ATGGATAAGCCACCGTAGTTATAACGTCCCCGCCCCCCACGTCGCATTTACAACCGC','ATGGATAAGTAGTAATATTTACCACGTCCCCGCCCCCCACGTCGCATTTACAACCGC','ATGGATAAGTGACCGTAATTACCACGTCCCCGCCCCCCACGTCGCATTTACAACCGC','ATGGATAAGCCACCGTAGTAACCACGTCCCCGCCCCCCACGTCGCATTTACAACCGC','ATGGATAAGCCATAGTATTAACCACGTCCCCGCCCCCCACGTCGCATTTACAACCGC','ATGGATAAGCCACCGTATTAATAACGTCCCCGCCCCCCACGTCGCATTTACAACCGC','ATGGATAAGCCACCGTATTAGCCATGACCCCGCCCCCCACGTCGCATTTACAACCGC','ATGGATAAGCCACCGTATTTATGATAACCCCGCCCCCCACGTCGCATTTACAACCGC']
    gp0=12.5
    gp_j={}
    for key in a:
        if key not in dic.keys():
            pass
        else:
            gp_j[key]=np.zeros(15)
            for i in range(15):
                gp_j[key][i]=gp0-math.log(dic[key][i],2) 
    return gp_j

g_p_dic = g_p(normalized_dic)  #bulk number of doubling times for innactive variants

def gp_lst():
    lst=[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]
    for i in range(15):
        sum_gp=0
        count=0
        for key in g_p_dic.keys():
            sum_gp+=g_p_dic[key][i]
            count+=1
        lst[i]=sum_gp/count
    return(lst)

#g_p_lst = gp_lst()     #bulk number of doubling times for each condition

g_p_lst = [12.580088912744081,
 12.5,
 11.62546921184205,
 10.859832270666054,
 10.42569687864241,
 11.916575412776796,
 12.5,
 11.495134752814488,
 10.989455130187363,
 10.440441756536798,
 12.485865789607656,
 12.5,
 12.130804792506707,
 11.266640146326294,
 10.813554599770724]  #g_p_lst is from amino acid slope analysis (see G1_AA_slope_from_read_counts.py) because it holistically captures Gp by combining synonymous negative control read counts.

def growth(): 
    new_dic={}
    t=430
    for key in normalized_dic.keys():
        new_dic[key]=np.zeros(15)
        for i in range(15):
                new_dic[key][i]=1/(t*math.log(math.e,2))*(math.log(normalized_dic[key][i],2)+g_p_lst[i])
    return(new_dic)
    
growth_rt_dic = growth() #growth rate for each variant at each induction condition and replicate

def potencymetric(dic):
    new_dic={}
    for key in dic.keys():
        new_dic[key]=np.zeros(6)
        new_dic[key][0]=np.polyfit([0,0.15,0.30,0.50],dic[key][1:5],1)[0]    #slope rep. 1
        new_dic[key][1]=np.polyfit([0,0.15,0.30,0.50],dic[key][6:10],1)[0]   #slope rep. 2
        new_dic[key][2]=np.polyfit([0,0.15,0.30,0.50],dic[key][11:],1)[0]    #slope rep. 3
        new_dic[key][3]=np.polyfit([0,0.15,0.30,0.50],dic[key][1:5],1)[1]    #y intercept rep. 1
        new_dic[key][4]=np.polyfit([0,0.15,0.30,0.50],dic[key][6:10],1)[1]   #y intercept rep. 2
        new_dic[key][5]=np.polyfit([0,0.15,0.30,0.50],dic[key][11:],1)[1]    #y intercept rep. 3
    return(new_dic)

slope_dic = potencymetric(growth_rt_dic) #slope and y-intercept dictionary

def average_metric(dic):
    ave_dic={}
    for key in dic.keys():
        ave_dic[key]=np.zeros(2)
        ave_dic[key][0]=(dic[key][0]+dic[key][1]+dic[key][2])/3    #slope
        ave_dic[key][1]=(dic[key][3]+dic[key][4]+dic[key][5])/3    #y intercept
    return(ave_dic)
    
def std_metric(dic):
    std_dic={}
    for key in dic.keys():
        std_dic[key]=np.zeros(2)
        std_dic[key][0]=stat.stdev([dic[key][0],dic[key][1],dic[key][2]])  #slope std dev
        std_dic[key][1]=stat.stdev([dic[key][3],dic[key][4],dic[key][5]])  #y intercept std dev
    return(std_dic)

average_dic = average_metric(slope_dic)  #dictionary of average slope and y-intercept across three repcliates

std_dic = std_metric(slope_dic) #dictionary of slope and y-intercept standard deviation across three repcliates

def significance_test(av_dic,st_dic):
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
    for key in average_dic.keys():
        if F_test(st_dic[key][0],st_dic['ATGGATAAGCCACCGTATTTACCACGTCCCCGCCCCCCACGTCGCATTTACAACCGC'][0]):  
            #t-test with unequal variance
            t=(av_dic[key][0]-av_dic['ATGGATAAGCCACCGTATTTACCACGTCCCCGCCCCCCACGTCGCATTTACAACCGC'][0])/(st_dic[key][0]**2/3+st_dic['ATGGATAAGCCACCGTATTTACCACGTCCCCGCCCCCCACGTCGCATTTACAACCGC'][0]**2/3)**0.5
            dof=(st_dic[key][0]**2/3+st_dic['ATGGATAAGCCACCGTATTTACCACGTCCCCGCCCCCCACGTCGCATTTACAACCGC'][0]**2/3)**2/((st_dic[key][0]**2/3)**2/2+(st_dic['ATGGATAAGCCACCGTATTTACCACGTCCCCGCCCCCCACGTCGCATTTACAACCGC'][0]**2/3)**2/2)
            p_value=tstat.cdf(t,dof)
            new_dic[key]=np.array([p_value,av_dic[key][0],av_dic[key][1],st_dic[key][0],st_dic[key][1],1])  
        else:
            #t-test with equal variance
            new_dic[key]=np.zeros(3)
            s_pooled=((2*st_dic[key][0]**2+2*st_dic['ATGGATAAGCCACCGTATTTACCACGTCCCCGCCCCCCACGTCGCATTTACAACCGC'][0]**2)/4)**0.5
            t=(av_dic[key][0]-av_dic['ATGGATAAGCCACCGTATTTACCACGTCCCCGCCCCCCACGTCGCATTTACAACCGC'][0])/s_pooled*(9/6)**0.5   
            dof=6-2
            p_value=tstat.cdf(t,dof)
            new_dic[key]=np.array([p_value,av_dic[key][0],av_dic[key][1],st_dic[key][0],st_dic[key][1],0])
    return(new_dic)

t_values = significance_test(average_dic,std_dic)   #1-tailed t-test relative to parental 

def p_value_array(dic):
    sort_lst=[]
    for i in dic.keys():
        sort_lst+=[(i,dic[i][0],dic[i][1],dic[i][2],dic[i][3],dic[i][4])]  #[name,p_value,aveslope,aveyint,stdevslope,stdevyint]
    return(sort_lst)

p_unsorted = p_value_array(t_values)  #list of p values

def sort(values):
    dtype=[('name', 'U60'), ('pvalue', 'U30'),('average_slope', 'U30'),('average_yint', 'U30'),('std_dev_slope', 'U30'),('std_dev_yint', 'U30')]   #originally said ('std_dev_yint', float)
    a = np.array(values, dtype=dtype)  
    return(np.sort(a, order='pvalue'))

sorted_values = sort(p_unsorted) #sort mutants by 1-tailed p-value

nameHandle = open('Gen_1_results_for_DNA.csv', 'w')  #writing .csv file
string = 'Sequence,p-value (1-tailed),Average slope (mM^-1min^-1),Average y-intercept (min^-1),Slope standard deviation (mM^-1min^-1),Y-intercept standard deviation (min^-1)'
nameHandle.write(string + '\n')
for i in range(len(sorted_values)):
    row = sorted_values[i]
    name = row[0]
    pval = row[1]
    slope = row[2]
    yint = row[3]
    slope_stdev = row[4]
    yint_stdev = row[5]
    string = name + ',' + pval  + ',' + slope  + ',' + yint  + ',' + slope_stdev  + ',' + yint_stdev
    nameHandle.write(string + '\n')
nameHandle.close()
