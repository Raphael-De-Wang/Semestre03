#!env python

import numpy as np
import matplotlib.pyplot as plt

def sample_call_rate(chr_file_name) :
    info_indi_list = []
    genotype_list  = []
    with open(chr_file_name) as handler :
        for line in handler :
            line = line.strip().split()
            info_indi = line[:6]
            genotype  = line[6:]
            info_indi_list.append(info_indi)
            genotype_list.append(genotype)
            
    allels1_ind    = np.arange(len(genotype)/2)*2
    allels2_ind    = np.arange(len(genotype)/2)*2+1    
    call_rate_list = np.sum((np.array(genotype_list)[:,allels1_ind]=='0') & (np.array(genotype_list)[:,allels2_ind]=='0'),axis=1) * 1. / len(allels1_ind)
    return (info_indi_list, genotype_list, call_rate_list)

def call_rate_hist(call_rate_list, bins=20):
    plt.hist(call_rate_list, bins=bins)
    plt.show()


    
(info_indi_list, genotype_list, call_rate_list) = sample_call_rate("../Data/chr2.ped")
call_rate_hist(call_rate_list)

(info_indi_list, genotype_list, call_rate_list) = sample_call_rate("../Data/chr5.ped")
call_rate_hist(call_rate_list)

(info_indi_list, genotype_list, call_rate_list) = sample_call_rate("../Data/chr13.ped")
call_rate_hist(call_rate_list)

(info_indi_list, genotype_list, call_rate_list) = sample_call_rate("../Data/chr16.ped")
call_rate_hist(call_rate_list) 

