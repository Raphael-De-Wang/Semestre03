#!env python

import sys
import argparse
import numpy as np


class hit(object) :
    def __init__(self,tab_line):
        # Fields: query id, subject id, % identity, alignment length, mismatches,
        #         gap opens, q. start, q. end, s. start, s. end, evalue, bit score
        (self.qid,self.sid,self.iden_percent, self.align_len, self.mismat,
         self.gaps, self.qstart, self.qend, self.sstart, self.send,
         self.evalue, self.bit_score) = tab_line.split('\t')
        self.iden_percent = float(self.iden_percent)
        self.align_len    = int(self.align_len)
        self.mismat       = int(self.mismat)
        self.gaps         = int(self.gaps)
        self.qstart       = int(self.qstart)
        self.qend         = int(self.qend)
        self.sstart       = int(self.sstart)
        self.send         = int(self.send)
        self.evalue       = float(self.evalue)
        self.bit_score    = float(self.bit_score)
        self.tab_line     = tab_line

    def is_combination_possible(self,contigs_lens_dict):
        qlen = contigs_lens_dict[self.qid]
        slen = contigs_lens_dict[self.sid]
        def is_terminus(spos,epos,seqlen):
            if spos == 1 or spos == seqlen or epos == 1 or epos == seqlen :
                return True
            else:
                return False
        if self.qid == self.sid :
            return False
        elif is_terminus(self.qstart,self.qend,qlen) and is_terminus(self.sstart,self.send,slen) :
            return True
        else:
            return False


def interface_standard():
    parser = argparse.ArgumentParser("Simple command for select best blastp hits.")
    parser.add_argument('-f','--hits-file-name', dest='blasthits', help='input blastp hits file name',required=True)
    parser.add_argument('-o','--output-file-name', dest='ofname', help='output blastp hits file name',required=True)
    parser.add_argument('-i','--iden-percent', dest='iden_percent', help='% identity')
    parser.add_argument('-a','--align-len', dest='align_len', help='alignment length')
    parser.add_argument('-m','--mismat', dest='mismat', help='mismatches')
    parser.add_argument('-g','--gaps', dest='gaps', help='gap opens')
    parser.add_argument('-e','--evalue', dest='evalue', help='evalue')
    parser.add_argument('-b','--bit-score', dest='bit_score', help='bit score')
    return parser.parse_args()

        
def load_blastn_tab_fname_to_obj_list(blastn_tab_fname):
    hits_list = []
    handler = open(blastn_tab_fname)
    for line in handler:
        if len(line) == 0 or line[0] == '#' :
            continue
        h = hit(line.strip())
        hits_list.append(h)
    return hits_list


def args2filter(args):
    if args.iden_percent :
        iden_percent = float(args.iden_percent)
    else :
        iden_percent = 0.
        
    if args.align_len :
        align_len = int(args.align_len)
    else :
        align_len=0
        
    if args.mismat :
        mismat = int(args.mismat)
    else :
        mismat = np.inf
        
    if args.gaps :
        gaps = int(args.gaps)
    else :
        gaps = np.inf

    if args.evalue :
        evalue = float(args.evalue)
    else :
        evalue = np.inf

    if args.bit_score :
        bit_score = float(args.bit_score)
    else :
        bit_score = np.inf
        
    return (iden_percent, align_len, mismat, gaps, evalue, bit_score)

def best_hits_filter(hits_list, iden_percent = 0., align_len=0, mismat=np.inf, gaps=np.inf, evalue=np.inf, bit_score=np.inf) :
    print iden_percent, align_len, mismat, gaps, evalue, bit_score
    best_hits_list = []
    for h in hits_list :
        if h.iden_percent > iden_percent and h.align_len > align_len and h.mismat < mismat and h.gaps < gaps and h.evalue < evalue and h.bit_score < bit_score :
            best_hits_list.append(h)
    return best_hits_list

def persist_best_hits(ofname, best_hits_list):
    handler = open(ofname,'w')
    for h in best_hits_list :
        handler.write(h.tab_line+'\n')
    handler.close()

if __name__ == "__main__":
    args = interface_standard()
    hits_list = load_blastn_tab_fname_to_obj_list(args.blasthits)
    best_hits_list = best_hits_filter(hits_list,*args2filter(args))
    persist_best_hits(args.ofname, best_hits_list)
