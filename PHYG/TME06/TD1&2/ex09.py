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
    parser.add_argument('-if', '--input-file-name', dest='ifname', help='input analysis name',required=True)
    parser.add_argument('-p', '--plot-file-name', dest='pname', help='plot file name',required=True)
    return parser.parse_args()


def load_fasta_file(fname):
    seqRec_list = []
    for seq_record in SeqIO.parse(fname, "fasta"):
        seqRec_list.append(seq_record)
    return seqRec_list

def plot()
    pass
        
if __name__ == "__main__":
    
    args = interface_standard()
    
    # ????
        
    plot()
        
