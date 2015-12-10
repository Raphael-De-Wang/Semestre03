#!env python

import csv
import argparse
import numpy as np
import matplotlib.pyplot as plt

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_protein


class cluster(object) :
    def __init__(self, cluster_name, record_lines):
        self.records = record_lines
        self.cluster_name = cluster_name
        self.gene_list = []
        self.fam_count = {}
        for record in record_lines :
            gene = record.split()[2]
            gene = gene.replace(">","")
            gene = gene.replace(".","")
            self.gene_list.append(gene.strip())
            fam  = gene.split("_")[-1]
            if self.fam_count.has_key(fam) :
               self.fam_count[fam] += 1
            else :
               self.fam_count[fam] = 1 


def interface_standard():
    parser = argparse.ArgumentParser("")
    parser.add_argument('-c', '--cluster-file-name', dest='cfname', help='cluster file name',required=True)
    parser.add_argument('-o', '--output-file-name', dest='oname', help='output file name',required=True)
    return parser.parse_args()


def build_comm_total_table (tab, families, cluster_list) :
    for i in range(1,len(families)+1) :
        comm_set = set()
        total_set= set()
        self_set = set()
        for c in cluster_list :
            cflag = True
            for fam in families[:i] :
                if fam in c.fam_count.keys() :
                    total_set.add(c.cluster_name)
                else :
                    cflag = False
                    continue
            if cflag :
                comm_set.add(c.cluster_name)
            if families[i-1] in c.fam_count.keys() :
                self_set.add(c.cluster_name)
        print i,families[i-1],len(comm_set),len(total_set),len(self_set)
        tab.append([i,families[i-1],len(comm_set),len(total_set),len(self_set)])
    return tab


def write_csv(csv_fname,table):
    with open(csv_fname, 'wb') as csvfile:
        spamwriter = csv.writer(csvfile, delimiter=',')
        for line in table :
            spamwriter.writerow(line)


if __name__ == "__main__":
    def encap_cluster(cluster_name, cluster_records, cluster_list):
        if cluster_name and len(cluster_records) :
            cluster_list.append(cluster(cluster_name, cluster_records))
        
    args = interface_standard()
    
    table = [["id","family","common","total","self"]]
    cluster_list = []
    cluster_name = None
    cluster_records = []
    with open(args.cfname) as handler :
        for line in handler :
            if line[0] == ">" :
                encap_cluster(cluster_name, cluster_records, cluster_list)
                cluster_name = line.strip()[1:]
                cluster_records = []
                continue
            cluster_records.append(line)
        encap_cluster(cluster_name, cluster_records, cluster_list)

    fam_set = set()
    for c in cluster_list :
        for fam in c.fam_count.keys() :
            fam_set.add(fam)
            
    print len(cluster_list)
    build_comm_total_table (table, list(fam_set), cluster_list)

    write_csv(args.oname,table)
