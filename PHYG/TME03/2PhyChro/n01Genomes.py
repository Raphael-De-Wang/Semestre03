#!/usr/bin/env python
# -*- coding: utf-8 -*-
###############################################################################
# Some small functions to get the characteristics of the chromosomes of a genomes

import re, os

## RETURN 3 lists: the names, the centromere positions, the last feature numbers
# for a given species
def chromosomesLists(path, speciesName):
   chromoFile=open(path+speciesName+'.ch','r')
   line=chromoFile.readline()
   ligneSplit=re.split('\t',line)
   name=ligneSplit[:-1]
   line=chromoFile.readline()
   ligneSplit=re.split('\t',line)
   centro=map(int,ligneSplit[:-1])
   line=chromoFile.readline()
   ligneSplit=re.split('\t',line)
   end=map(int,ligneSplit[:-1])
   chromoFile.close
   return name,centro,end

# Return the list of the different .def found in the 01Genomes/ directory
def ListShortName(path):
   ## Genome
   genomesL=os.listdir(path)
   genomesL.sort()
   shortNameList=[]
   for f in genomesL:
      if re.match('(.*)\.def$',f):
         res=re.match('(.*)\.def$',f)
         shortNameList.append(res.group(1))
   return shortNameList



