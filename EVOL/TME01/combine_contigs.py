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
    parser.add_argument('-o', '--order', metavar='N', type=int, dest='contigs_order', nargs='+', help="contigs order, N relate to an integer, positive or negative relate to the direction of strain. ",required=True)
    parser.add_argument('-s', '--leading-lagging', metavar='S', type=int, dest='strains', nargs='+', help='leading(0) or lagging(other) strain', required=True)
    parser.add_argument('-L','--overlap', metavar='L', type=int, dest='overlaps', nargs='+', help="length of overlaps between contigs",required=True)
    parser.add_argument('-n','--combine-name', dest='cname',help="Give a name to combined contigs")
    return parser.parse_args()
    
def load_input_fasta_to_list(fname):
    seqRec_list = []
    for seq_record in SeqIO.parse(fname, "fasta"):
        seqRec_list.append(seq_record)
    return seqRec_list

def save_combined_contigs_to_fasta(fname,seq_record):
    SeqIO.write([seq_record], fname, "fasta")

def get_arranged_contig_by_order(seqRec_list, order, leading):
    seq_record = seqRec_list[abs(order)-1]
    if order < 0 :
        seq_record = seq_record[::-1]
    if not leading :
        seq_record = seq_record.reverse_complement()
    return seq_record
    
def combine_contigs(seqRec_list,order_list,overlap_list,strain_list,cname=None):
    seq = ""
    overlap_seq = ""
    ol_list = [0] + list(overlap_list) + [0] # to insert a zero, there is no overlap on first contigs
    if cname == None:
        cname = "gene"
    for ind, order in enumerate(order_list) : 
        overlap = int(ol_list[ind])
        seq_record = get_arranged_contig_by_order(seqRec_list, order, strain_list[ind]==0)
        if cmp(overlap_seq,seq_record.seq[:overlap].tostring()) <> 0 :
            print overlap_seq
            print seq_record.seq[:overlap]
            raise ValueError("Overlap in contigs meet conflicts, contig id [%s]"%order)
        seq += seq_record.seq[overlap:]
        overlap = ol_list[ind+1]
        overlap_seq = seq_record.seq[-overlap:].tostring()
    return SeqRecord(seq,id=cname,name=cname,description="|combine contigs|")

args = interface_standard()
seqRec_list = load_input_fasta_to_list(args.ifname)
seq_record = combine_contigs(seqRec_list,args.contigs_order,args.overlaps,args.strains,args.cname)
save_combined_contigs_to_fasta(args.ofname,seq_record)

