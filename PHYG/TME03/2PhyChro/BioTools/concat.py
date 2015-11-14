#!/usr/bin/python
# -*- coding: Latin-1 -*-

import sys
import os

# Ingrid Lafontaine - mars 2006 (last update Cecile Neuveglise 30 juin 2008)

###########################################
# usage : ./concat.py fichier_fasta

# fichier_fasta contient toutes les sequences au format fasta multiple
###########################################

    
#------------------------------------------ 
def Fastaout(nam,alin,nb) :
    sortie = nam+".fasta"
    print sortie
    filo = open(sortie,'w')
    filo.write(" "+str(len(alin.keys()))+" "+str(len(alin[alin.keys()[0]])-nb)+"\n")
    for k in alin.keys():
        filo.write(">"+k+"\t"+alin[k])
    filo.close()

#-----------------------------------------
def main():
    file=open(sys.argv[1],"r")
    seq={}

    line = file.readline()
    nseq = ''
    nbN=0
    while line :
        if line[0] == ">":
            nseq = line.split()[0][1:]
	    if nseq not in seq.keys():
		seq[nseq]=""
		print nseq
        elif line[:-1] != '': ## empty line
	    for l in line.split():                                         
            	seq[nseq] += l
	    seq[nseq] += '\n'
	    nbN+=1
        line = file.readline()
    print nbN,len(seq.keys()),nbN/len(seq.keys())
    Fastaout(sys.argv[1],seq,nbN/len(seq.keys()))
    file.close()

if len(sys.argv) < 1 :
    print "\n usage : ./concat.py fichier_fasta1\n"
else :
    main()                                                                                            
