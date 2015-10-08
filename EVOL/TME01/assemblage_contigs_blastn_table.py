#!env python

import sys
import numpy as np

class align_record(object) :
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

def load_blastn_tab_fname_to_obj_list(blastn_tab_fname):
    record_list = []
    handler = open(blastn_tab_fname)
    for line in handler:
        if len(line) == 0 or line[0] == '#' :
            continue
        ar = align_record(line.strip())
        record_list.append(ar)
    return record_list

def exact_contigs_len_list(blastn_record_list):
    contigs_lens_dict = {}
    for ar in blastn_record_list :
        if ar.qid == ar.sid and ar.iden_percent >= 1. :
            if not contigs_lens_dict.has_key(ar.qid) or contigs_lens_dict[ar.qid] < ar.align_len :
                contigs_lens_dict[ar.qid] = ar.align_len
    return contigs_lens_dict

def exact_assemblage_possible_list(blastn_record_list, contigs_lens_dict):
    ass_pair_list = []
    for ar in blastn_record_list :
        if ar.is_combination_possible(contigs_lens_dict):
            if not (ar.qid,ar.sid) in ass_pair_list and not (ar.sid,ar.qid) in ass_pair_list :
                ass_pair_list.append((ar.qid,ar.sid))
                print ar.qid, ar.qend > ar.qstart, ar.sid, ar.send > ar.sstart, abs(ar.send - ar.sstart) + 1
    return ass_pair_list

fname = sys.argv[1]
record_list       = load_blastn_tab_fname_to_obj_list(fname)
contigs_lens_dict = exact_contigs_len_list(record_list)
ass_pair_list     = exact_assemblage_possible_list(record_list, contigs_lens_dict)

