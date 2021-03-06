#!env python

import argparse
import numpy as np
import matplotlib.pyplot as plt

from matplotlib.backends.backend_pdf import PdfPages

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


def plot(prot_num_list, strain_name):
    fig, ax = plt.subplots()
    plt.hist(prot_num_list)
    plt.title(strain_name)
    plt.xlabel("protein size")
    plt.ylabel("distribution")
    

        
if __name__ == "__main__":
    
    args = interface_standard()
    
    plot_handle = PdfPages(args.pname)
    
    for fname in args.ifname :
        prot_size_list = []
        seq_list = load_fasta_file(fname)
        for seq in seq_list :
            prot_size_list.append(len(seq))
        plot(prot_size_list, fname.split("/")[-1][:-14])
        plot_handle.savefig()
        
    plot_handle.close()

