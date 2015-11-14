#!/usr/bin/env python
# -*- coding: utf-8 -*-
## Script to convert muscle-gb files into phylip files

import re,os,sys
import n01Genomes


def _main_(argv):

   if len(argv)!=4:
      print "4 arguments are expected:\n\
      -> the clade Name (as 'Yeast')\n\
      -> delta\n\
      -> 0 (to built the whole tree) / 1 (to choose)\n"
      exit(1)
   clade=argv[1]
   delta=argv[2]
   path='../../'+clade+'/12OrthFamilies/Delta'+delta+'/'
   Qall=int(argv[3])
   
   ### va chercher les differents genomes
   iniShortNameList=n01Genomes.ListShortName('../../'+clade)
   iniShortNameList.sort()

   if Qall!=0:   ## Question/Interface
      i=0
      for x in iniShortNameList:
         i+=1
         print i,'\t',x

      questionGenome=raw_input("Which genomes do you want to align? ... write their numbers separated by spaces \n")
      questionLigne=re.split(" ",questionGenome)
      shortNameList=[]
      for num in questionLigne: #questionLigne[:-1]:
            if int(num)<=len(iniShortNameList):
               shortNameList.append(iniShortNameList[int(num)-1])
            else:
               print "!!! All number has to correspond to a genome  !!!!"
   else:
      shortNameList=iniShortNameList

   shortNameList.sort()
   nbg=len(shortNameList)
   print 'Alignement of the species', shortNameList, ' \n'

   output=open(path+''.join(map(lambda x:x[0],shortNameList))+'_AlignSeq.data','w')

   listFiles=os.listdir(path)
   seq={} ## le dico des sequences {'LAKL': 'ANGPIPDARN...',..}
   for shName in shortNameList:
      seq[shName]=''
   
   i=1
   si='%05d'%i
   while 'Prot_'+si+'.muscle-gb' in listFiles:
      fileI=open(path+'Prot_'+si+'.muscle-gb','r')
      for ligne in fileI:
         if re.match('>',ligne):
            lList=re.split('\t',ligne)
            shName=lList[0][1:5]
            if shName in shortNameList:
               interesting=1
            else:
               interesting=0
         else:
            if interesting:
               lList=re.split(' ',ligne)
               seq[shName]+=''.join(lList)[:-1]
      fileI.close()
      i+=1
      si='%05d'%i
   exactLseq=len(seq[shortNameList[0]])
   output.write('   '+str(nbg)+'   '+str(exactLseq)+'\n')
   seqC={}
   for shName in shortNameList:
      #print shName,seq[shName]
      seqC[shName]=shName+'        '
      for i in range(exactLseq/10):
         seqC[shName]+=seq[shName][i*10:10*(i+1)]+' '
         if (i+1)%6==0:
            seqC[shName]+='\n            '
      seqC[shName]+=seq[shName][(i+1)*10:]+' \n'
   for i in range(len(seqC[shortNameList[0]])/79):
      for shName in shortNameList:
         output.write(seqC[shName][i*79:79*(i+1)])
      output.write('\n')
   for shName in shortNameList:
      output.write(seqC[shName][(i+1)*79:])
   output.close()

#_main_([0,'Yeast','6','0'])
#_main_([0,'Vertebrate','6','0'])
if __name__=="__main__":
   _main_(sys.argv)
