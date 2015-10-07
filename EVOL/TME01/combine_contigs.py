#!env python

import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_protein

def interface_standard():
    parser = argparse.ArgumentParser("Simple command for combining specified contigs")
    parser.add_argument('-in', '--input-file-name', dest='ifname', help='input fasta name',required=True)
    parser.add_argument('-on', '--output-file-name', dest='ofname', help='output fasta name',required=True)
    parser.add_argument('-o', '--order', metavar='N', dest='contigs_order', nargs='+', help="contigs order, N relate to an integer, positive and negative relate to the direction of strain. With ` means laggings strain. Ex. : 3 means contig3 plus direction 5`-3`, 5` means contig5 minus direction 5`-3`.",required=True)
    parser.add_argument('-l','--overlap', metavar='L', dest='overlaps', nargs='+', help="length of overlaps between contigs",required=True)
    parser.add_argument('-n','--combine-name', dest='cname',help="Give a name to combined contigs")
    return parser.parse_args()
    
def load_input_fasta_to_list(fname):
    seqRec_list = []
    for seq_record in SeqIO.parse(fname, "fasta"):
        seqRec_list.append(seq_record)
    return seqRec_list

def save_combined_contigs_to_fasta(fname,seq_record):
    SeqIO.write([seq_record], fname, "fasta")

def combine_contigs(seqRec_list,order_list,overlap_list,cname=None):
    seq = ""
    overlap_seq = ""
    if cname == None:
        cname = "gene"
    for seq_record, order, overlap in zip(seqRec_list, order_list, [0]+list(overlap_list)) : # to insert a zero, there is no overlap on first contigs
        if cmp(overlap_seq,seq_record.seq[:overlap]) <> 0 :
            raise ValueError("Overlap in contigs meet conflicts, contig id [%s]"%order)
        seq += seq_record.seq[overlap:]
    return SeqRecord(seq,id=cname,name=cname,description="combine contigs")

args = interface_standard()
seqRec_list = load_input_fasta(args.ifname)
seq_record = combine_contigs(seqRec_list,args.contigs_order,args.overlaps,args.cname)
save_combined_contigs_to_fasta(fname,seq_record)