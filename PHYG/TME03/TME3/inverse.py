#!env python

import numpy as np

centremers_list = [(301,302), (1285,1286), (1711,1712), (2524,2525), (3557,3558),
                (4475,4476),  (5301,5302), (5849,5850), (6430,6431), (7745,7746),
            (8005,8006), (8978,8979), (9938,9939), (10401,10402), (11535,11536),
            (12399,12400), (13265,13266), (13366,13367), (14340,14341), (15224,15225),
            (16031,16032), (16594,16595), (17646,17647)]

def read_fake_fasta(fname, cen_pos):
    handler = open(fname,'r')
    chrs_list = []
    for line in handler :
        if line[0] == '>':
            gname = line[1:]
            continue
        chrs_list.append([ int(c) for c in line.split()[:-1] ])
    return (gname, insert_centremers(cen_pos,chrs_list))

def insert_0_chr(centremer_pos, chromosome):
    for i in range(len(chromosome)):
        if chromosome[i:i+2] == list(centremer_pos) :
            chromosome.insert(i,0)
            return chromosome
    raise ValueError("Centremer Position unfound [%d,%d]"%centremer_pos)

def insert_centremers(centremers_list,chromosomes_list):
    return [ insert_0_chr(centremers_list[i], chm) for i,chm in enumerate(chromosomes_list) ]
        
def inversion(centremer,chromosome):
    def not_cross_cntm():
        pass
    def rand_pos():
        rand_chr  # uniform
    def rand_num_gen():
        rand_gen = np.random.poisson(5) # poisson
    def rand_inverse():
        pass
    return rand_inverse(cntm_pos,chm)

def fission(centremers_list,chromosomes_list):
    def rand_chr(chrs_list):
        return int(np.floor(np.random.uniform(len(chrs_list))))
    
    def rand_gene_pos(chromosome):
        return int(np.floor(np.random.uniform(len(chromosome))))

    def chr_breaker(chromosome):
        brk_pos = rand_gene_pos(chromosome)
        return (chromosome[:brk_pos],chromosome[brk_pos+1:])
    
    def has_cntm(brk_chr):
        return sum(np.array(brk_chr) == 0) > 0
    
    def rand_insert_cntm(brk_chr):
        pos = int(np.floor(np.random.uniform(1,len(brk_chr))))
        brk_chr.insert(pos, 0)
        return brk_chr
    
    def update_chr(chrs_list, brk_chr1, brk_chr2):
        if not has_cntm(brk_chr1) :
            brk_chr1 = rand_insert_cntm(brk_chr1)
        elif not has_cntm(brk_chr2) :
            brk_chr2 = rand_insert_cntm(brk_chr2)
        chrs_list.append(brk_chr1)
        chrs_list.append(brk_chr2)
        return chrs_list
    
    chr_num = rand_chr(chromosomes_list)
    chromo  = chromosomes_list.pop(chr_num)
    (chr1, chr2) = chr_breaker(chromo)
    return update_chr(chromosomes_list, chr1, chr2)
    
(gname, chrs) = read_fake_fasta("ancestor.txt", centremers_list)

print len(chrs)

chl = fission(centremers_list,chrs)

print len(chl)

print sum([ np.sum(np.array(c)==0) for c in chl ])


