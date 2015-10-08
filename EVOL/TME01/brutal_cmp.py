#!env python

import argparse

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_protein

def interface_standard():
    parser = argparse.ArgumentParser("Simple command for compaire ends of contigs. ")
    parser.add_argument('-in', '--fasta-file-name', dest='ifname', help='input fasta name',required=True)
    parser.add_argument('-c1','--contig1-name', dest='contig1',help='contig 1 name',required=True)
    parser.add_argument('-c2','--contig2-name', dest='contig2',help='contig 2 name',required=True)
    parser.add_argument('-s','--seuil', type=int, dest='seuil',help='max number of nu compared from each contig end. Default: 50')
    return parser.parse_args()

def load_input_fasta_to_dict(fname):
    seqRec_dict = {}
    for seq_record in SeqIO.parse(fname, "fasta"):
        seqRec_dict[seq_record.id] = seq_record
    return seqRec_dict

def iterative_end_pair(head,tail):
    if cmp(head,tail) == 0 :
        return len(head)
    for i in range(2,len(head)) :
        if cmp(head[:i],tail[-i:]) == 0 :
            print tail[-i:]
            return len(tail[-i:])
    return -1
    
def compaire_4_pair(contig1,contig2,seuil=None):
    if seuil is None :
        seuil = 50
    i = iterative_end_pair(contig1[:seuil].seq.tostring(),contig2[-seuil:].seq.tostring())
    if i >= 0 :
        return (contig1.id,contig2.id,i)
    i = iterative_end_pair(contig2[:seuil].seq.tostring(),contig1[-seuil:].seq.tostring())
    if i >= 0 :
        return (contig2.id,contig1.id,i)
    i = iterative_end_pair(contig1.reverse_complement()[:seuil].seq.tostring(),contig2[-seuil:].seq.tostring())
    if i >= 0 :
        return (contig1.id+" compl",contig2.id,i)
    i = iterative_end_pair(contig2[:seuil].seq.tostring(),contig1.reverse_complement()[-seuil:].seq.tostring())
    if i >= 0 :
        return (contig2.id,contig1.id+" compl",i)
    
args = interface_standard()
seqRec_dict = load_input_fasta_to_dict(args.ifname)
print compaire_4_pair(seqRec_dict[args.contig1],seqRec_dict[args.contig2],args.seuil)
