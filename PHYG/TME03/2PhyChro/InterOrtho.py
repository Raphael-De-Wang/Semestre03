#!/usr/bin/env python
# -*- coding: utf-8 -*-
## n01SyntheseBlocs

############################################################################
## Construction de GENE GENEC
## Construction de BLOCK avec signe, inclusion et overlap
## et deletion des blocks inclus n'apportant aucune ancre de plus ds les 2 genomes
############################################################################

import math,re,os,sys
import gb2phylip


#################################
##  MAIN      MAIN      MAIN   ##
# fonctions appelees du script2 #
#################################

def ListShortName(organismName):## au lieu de dicoListShortName
   ## Genome
   genomesL=os.listdir('../../'+organismName+'/01Genomes/')
   genomesL.sort()
   shortNameList=[]
   for f in genomesL:
      if re.match('([A-z]+)\.def',f):
         res=re.match('([A-z]+)\.def',f)
         shortNameList.append(res.group(1))
   return shortNameList
   #return ['ERGO','KLLA','LAKL']

######
## On lit '.prt'
def synthesePrt(name,group):
   fileO=open('../../'+group+'/01Genomes/'+name+'.prt','r')
   dicoGene={}#numFeat:(name,sequence)
   for ligne in fileO:
      if re.match('>',ligne):
         ligneList=re.split("\t",ligne)
         #print ligneList
         dicoGene[int(ligneList[8])]=[ligneList[0],'']
      else:
         dicoGene[int(ligneList[8])][1]+=ligne
      #print dicoGene[ligneList[9]][0],ligneList[9]
   fileO.close()
   return dicoGene

######
## On lit '.orth.pairs'
def synthesePairs(name1,name2,group,delta):
   dico1={}
   dico2={}
   fileO=open('../../'+group+'/11Blocks/Delta'+delta+
            '/OrthBlocks/'+name1+'.'+name2+'.orth.pairs','r')

   for ligne in fileO:
      ligneList=re.split(" ",ligne)
      num1=int(ligneList[4])
      num2=int(ligneList[9])
      percent=int(ligneList[10])
      if num1 in dico1.keys():
         dico1[num1].append((num2,percent))
      else:
         dico1[num1]=[(num2,percent)]
      if num2 in dico2.keys():
         dico2[num2].append((num1,percent))
      else:
         dico2[num2]=[(num1,percent)]

   fileO.close()
   return dico1,dico2

######
## On supprime si homologue non unique
def unicHomo(dico1,dico2):
   for k,v in dico1.copy().iteritems():
      if len(v)>1 or len(dico2[v[0][0]])>1:
         del dico1[k]
      else:
         dico1[k]=v[0]
   for k,v in dico2.copy().iteritems():
      if len(v)>1 or v[0][0] not in dico1.keys():
         del dico2[k]
      else:
         dico2[k]=v[0]
   return dico1,dico2

#############
## on fait l'intersection vraiment strict
def intersection(listDico,nbG):
   listOrtho=listDico[0][1].keys()
   listOrtho=map(lambda x:[x],listOrtho)

   listPercent=[]
   for x in listOrtho:
      listPercent.append(0)

   for i in range(1,nbG):
      #dico=listDico[i][i+1]
      for j in range(len(listOrtho)):
         sameOrthos=listOrtho[j]
         if len(sameOrthos)==i:
            if sameOrthos[0] in listDico[0][i].keys():
               h=listDico[0][i][sameOrthos[0]][0]
               p=listDico[0][i][sameOrthos[0]][1]
               for k,o in enumerate(sameOrthos[1:]):
                  if o in listDico[k+1][i]:
                     if listDico[k+1][i][o][0]==h:
                        p+=listDico[k+1][i][o][1]
                     else:
                        p=0
                        break
                  else:
                     p=0
                     break
               if p:
                  listOrtho[j].append(h)
                  listPercent[j]+=p
               else:
                  listPercent[j]=0

   return listOrtho,listPercent


#############
## on fait l'intersection
def premNeg(liste):
   for i in range(len(liste)):
      if liste[i]<0:
         liste[i]=-1*liste[i]
         return i+1
   return 0

def compO(x,l):
   liste=[-1*x]
   for i in range(l-1):
      liste.append(0)
   return liste

def intersectionPlusSouple(listDico,nbG):
   listOrtho=set(listDico[0][1].keys())
   for i in range(2,nbG):
      listOrtho=set.union(listOrtho,set(listDico[0][i].keys()))
   listOrtho=list(listOrtho)
   listOrtho=map(lambda x:compO(x,nbG),listOrtho)

   listPercent=[]
   for x in listOrtho:
      listPercent.append(0)
   nbComp=[]
   for x in listOrtho:
      nbComp.append(0)
      
   for j in range(len(listOrtho)):
      sameOrthos=listOrtho[j]
      g=premNeg(sameOrthos)
      while g:
         g-=1
         for i in range(1,nbG):
            if i!=g:
               if sameOrthos[g] in listDico[g][i].keys():
                  h=listDico[g][i][sameOrthos[g]][0]
                  p=listDico[g][i][sameOrthos[g]][1]
                  if sameOrthos[i]!=0:
                     if h!=abs(sameOrthos[i]):
                        listOrtho[j]=[]
                        break
                  else:
                     sameOrthos[i]=-1*h
         if listOrtho[j]==[]:
            break
         g=premNeg(sameOrthos)
   for j in range(len(listOrtho)):
      if len(listOrtho[j])==nbG and 0 not in listOrtho[j]:
         for g in range(0,nbG-1):
            for i in range(g+1,nbG):
               if listOrtho[j][g] in listDico[g][i].keys():
                  if listDico[g][i][listOrtho[j][g]][0]==listOrtho[j][i]:
                     #print listOrtho[k][i],listDico[i][j][listOrtho[k][i]][1]
                     listPercent[j]+=listDico[g][i][listOrtho[j][g]][1]
                     nbComp[j]+=1
                  else:
                     print "hhhhhh"
   return listOrtho,listPercent,nbComp

#############
## on ecrit ds les fichiers
def writeFiles(group,names,prtDico,listOrtho,listPercent,nbComp,miniP,delta):
   ## remove de previous prots
   protF=os.listdir('../../'+group+'/12OrthFamilies/Delta'+delta+'/')
   for i in range(len(protF)-1,-1,-1):
      if re.match('Prot_',protF[i]):
         os.remove('../../'+group+'/12OrthFamilies/Delta'+delta+'/'+protF[i])

   fileList=[]
   j=0
   for i in range(len(listPercent)):
      if nbComp[i]!=0:
         listPercent[i]=float(listPercent[i])/nbComp[i]/100
   while max(listPercent)>float(miniP)/100:
      maxiP=max(listPercent)
      index=listPercent.index(maxiP)
      listPercent[index]=-1
      ortho=listOrtho[index]
      if ortho!=[] and 0 not in ortho:
         f=open('../../'+group+'/12OrthFamilies/Delta'+delta+'/Prot_'+('%05d'% (j+1))+'.fasta','w')
         for i,x in enumerate(names):
            num=ortho[i]
            #print i,num
            f.write('>'+x+'   '+prtDico[i][num][0]+'  '+str(num)+'    '+str(maxiP)+'\n'+prtDico[i][num][1])
         j+=1
         f.close() 

# summary file
def writeFile(group,names,prtDico,listOrtho,listPercent,nbComp,miniP,delta):
   fileOut=open('../../'+group+'/12OrthFamilies/Delta'+delta+'/'+''.join(map(lambda x:x[0],names))+'_OrthoList.txt','w')
   fileOut.write('\t'.join(names)+'\tPercent'+miniP+'\n')
   fileList=[]
   j=0
   for i in range(len(listPercent)):
      if nbComp[i]!=0:
         listPercent[i]=float(listPercent[i])/nbComp[i]/100
   while max(listPercent)>float(miniP)/100:##for k in range(len(listPercent)):
      maxiP=max(listPercent)
      index=listPercent.index(maxiP)
      listPercent[index]=-1
      ortho=listOrtho[index]
      if ortho!=[] and 0 not in ortho:
         s='\t'.join(map(str,ortho))
         fileOut.write(s+'\t'+str(maxiP)+'\n')
         j+=1

   fileOut.close()



def _main_(argv):
   if len(argv)!=4:
      print "3 arguments are expected:\n -> the clade name (ex:'Yeast' or 'Vertebrate')\n -> a value of delta (in [1,6] for instance)\n -> the minimum average of similarity between orthologs of the same family\n"
      exit(1)
   group=argv[1]
   delta=argv[2]
   miniP=argv[3]
   os.system('mkdir ../../'+group+'/12OrthFamilies/')
   os.system('mkdir ../../'+group+'/12OrthFamilies/Delta'+delta+'/')
   nameList=ListShortName(group)
   #nameList=['CAAL','CADU','CAGL']
   nameList.sort()
   nbGenomes=len(nameList)
   print nbGenomes,nameList

   #### definir la list des dico des Genes [{num:(name,seq)},{},...]
   dicoGene=[]
   for name in nameList:
      dicoGene.append(synthesePrt(name,group))

   ####### definir la list de list des dico d'homo
   listlistHomo=[]
   for name in nameList:
      listlistHomo.append([])
   for i in range(nbGenomes-1):
      listlistHomo[i].append({})
      for j in range(i+1,nbGenomes):
         dico1,dico2=synthesePairs(nameList[i],nameList[j],group,delta)
         dico1,dico2=unicHomo(dico1,dico2)
         listlistHomo[i].append(dico1)
         listlistHomo[j].append(dico2)
   listlistHomo[nbGenomes-1].append({})

   ##nb de comparaison
   #som=nbGenomes*(nbGenomes-1)/2.*100.
   
   #listOrtho,listPercent=intersection(listlistHomo,len(nameList))
   #writeFiles(group,nameList,dicoGene,listOrtho,listPercent,som)
   listOrtho,listPercent,nbComp=intersectionPlusSouple(listlistHomo,len(nameList))
   writeFiles(group,nameList,dicoGene,listOrtho,listPercent[:],nbComp,miniP,delta)
   writeFile(group,nameList,dicoGene,listOrtho,listPercent[:],nbComp,miniP,delta)

   gb2phylip._main_([0,group,delta,0])

_main_(sys.argv)













