#!env python

import argparse
import numpy as np
import matplotlib.pyplot as plt

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_protein


def interface_standard():
    parser = argparse.ArgumentParser("")
    parser.add_argument('-if', '--input-file-name', dest='ifname', help='input fasta name',required=True, nargs='+')
    parser.add_argument('-p', '--plot-file-name', dest='pname', help='plot file name',required=True)
    return parser.parse_args()


def load_fasta_file(fname):
    seqRec_list = []
    for seq_record in SeqIO.parse(fname, "fasta"):
        seqRec_list.append(seq_record)
    return seqRec_list


def plot(prot_num_list, strain_name_list, plot_name, width=0.5) :
    fig, ax = plt.subplots(figsize=(12, 12))
    plt.subplots_adjust(bottom=0.35)
    ind = np.arange(len(prot_num_list))
    ax.bar(ind, prot_num_list, width, color='r')
    ax.set_title('the distribution of the number of predicted proteins by strain')
    ax.set_xticks(ind + width/2.)
    ax.set_xticklabels(strain_name_list,rotation=270)
    fig.savefig(plot_name)

        
if __name__ == "__main__":
    
    args = interface_standard()

    prot_num_list = []
    strain_name_list = []
    
    for fname in args.ifname :
        seq_list = load_fasta_file(fname)
        prot_num_list.append(len(seq_list))
        strain_name_list.append(fname.split("/")[-1][:-14])

    plot(prot_num_list, strain_name_list,args.pname)
        
        

