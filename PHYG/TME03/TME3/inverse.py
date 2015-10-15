#!env python

import numpy as np

"""
centremers_list = [(301,302), (1285,1286), (1711,1712), (2524,2526), (3557,3558),
                (4475,4476),  (5301,5302), (5849,5850), (6430,6431), (7745,7746),
            (8005,8006), (8978,8979), (9938,9939), (10401,10402), (11535,11536),
            (12399,12400), (13265,13266), (13366,13367), (14340,14341), (15224,15245),
            (16031,16032), (16594,16595), (17646,17647)]
"""

def insert_centremers(centremers_list,chromosomes_list):
    def insert_0_chr(centremer_pos, chromosome):
        return chromosome.insert(centremer_pos,0)
    return [ insert_0_chr(centremers_list[i], chm) for i,chm in enumerate(chromosomes_list) ]
        
def read_fake_fasta(fname, cen_pos):
    handler = open(fname)
    return insert_centremers(centremers_list,chromosomes_list)

def inversion(centremers_list,chromosomes_list):
    def not_cross_cntm():
        pass
    def gen_rand_chr():
        rand_chr  # uniform
    def gen_rand_gen_pos():
        rand_gen = np.random.poisson(5) # poisson
    def rand_inverse():
        pass
    return [rand_inverse(cntm_pos,chm) for cntm_pos,chm in zip(centremers_list,chromosomes_list)] 

def fission(centremers_list,chromosomes_list):
    def gen_rand_chr():
        rand_chr  # uniform
    def gen_rand_gen():
        rand_gen  # poisson
    def chr_breaker():
        pass
    def has_cntm():
        return True
    def rand_insert_cntm(): # warning : this change the numbers of chromosome, so centremers_list has to be updated as well 
        pass
    def update_chr():
        pass
    return (centremers_list+[new_pos],chromosomes_list+[new_chr])
    
def rand_evens
