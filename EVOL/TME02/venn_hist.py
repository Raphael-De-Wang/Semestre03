#!env python

import sys
import argparse
import pylab as plt
from matplotlib_venn import venn3, venn3_circles, venn2, venn2_circles

# Salmonella enterica -- CP
# Escherichia coli -- NC
# Staphylococcus aureus -- BX

def load_fam_silix(fname):
    fam_silix_dict = {}
    handler = open(fname)
    for line in handler :
        (fam,pid) = line.split()
        pid = pid.split('|')[1]
        if fam_silix_dict.has_key(fam) :
            fam_silix_dict[fam].append(pid)
        else :
            fam_silix_dict[fam] = [pid]
    handler.close()
    return fam_silix_dict

def venn_statis(fam_silix_dict):
    cp = 0
    nc = 0
    bx = 0
    cp_nc = 0
    cp_bx = 0
    nc_bx = 0
    cnb   = 0
    for fam,pids in fam_silix_dict.iteritems() :
        comp = "".join(pids)
        if "CP" in comp and "NC" in comp and "BX" in comp :
            cnb += 1
        elif "CP" in comp and "NC" in comp :
            cp_nc += 1
        elif "CP" in comp and "BX" in comp :
            cp_bx += 1
        elif "NC" in comp and "BX" in comp :
            nc_bx += 1
        elif "CP" in comp :
            cp += 1
        elif "BX" in comp :
            bx += 1
        elif "NC" in comp :
            nc += 1
        else :
            print pids
            raise ValueError("Missing Match family [%s]"%fam)
    return (cp,nc,bx,cp_nc,cp_bx,nc_bx,cnb)

def plot_venn3(cp,nc,bx,cp_nc,cp_bx,nc_bx,cnb,plotName=None) :
    fig = plt.figure()
    venn3(subsets={'100':cp,'010':nc,'001':bx,'110':cp_nc,'101':cp_bx,'011':nc_bx,'111':cnb},set_labels = ('Salmonella enterica', 'Escherichia coli', 'Staphylococcus aureus'))
    plt.title("Venn diagram")
    if plotName :
        plt.savefig(plotName)
    else :
        plt.show()

def interface_standard() :
    parser = argparse.ArgumentParser("Simple command for plot venns diagram.")
    parser.add_argument('-f','--fam-silix-file-name', dest='ifname', help='input family file name',required=True)
    parser.add_argument('-o','--plot-file-name', dest='ofname', help='output venn diagram plot file name')
    return parser.parse_args()
        
if __name__ == "__main__":
    args = interface_standard()
    fam_silix_dict = load_fam_silix(args.ifname)
    (cp,nc,bx,cp_nc,cp_bx,nc_bx,cnb) = venn_statis(fam_silix_dict)
    plot_venn3(cp,nc,bx,cp_nc,cp_bx,nc_bx,cnb,args.ofname) 
