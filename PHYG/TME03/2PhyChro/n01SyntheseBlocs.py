#!/usr/bin/env python
# -*- coding: utf-8 -*-
## n01SyntheseBlocs

############################################################################
## Construction de GENE GENEC
## Construction de BLOCK avec signe, inclusion et overlap
## et deletion des blocks inclus n'apportant aucune ancre de plus ds les 2 genomes
############################################################################

import re
from n00Objets import *

#################################
##  MAIN      MAIN      MAIN   ##
# fonctions appelees du script2 #
#################################

######
## On complete l'objet GENOME avec les features found in '.def'
def syntheseDef(genome,path):
   fileO=open(path+genome.name+'.def','r')
   fileO.readline()
   listFeature=[]
   for ligne in fileO:
      ligneList=re.split("\t",ligne)
      if ligneList[5]=='+':
         strand=1
      else:
         strand=-1
      if ligneList[6]=='t':
         sens=1
      else:
         sens=-1
      listFeature.append(Gene(ligneList[0],ligneList[1],genome,int(ligneList[2]),
                        (ligneList[3],ligneList[4]),strand,
                        sens,int(ligneList[9]),int(ligneList[8])))
   genome.setFeatures(listFeature)  
   fileO.close()

######
## On complete l'objet GENE en lui associant un GENEC avec des homologues BDBH ou non
# a partir du fichier '.orth.pairs'
def synthesePairs(genome1,genome2,path,i=0):
   if genome1.name<genome2.name:
      g1=genome1
      g2=genome2
   else:
      g2=genome1
      g1=genome2
   shName1=g1.name
   shName2=g2.name
   listFeatures1=g1.features
   listFeatures2=g2.features
     
   for g in listFeatures1:
      if g.genre=='gene':
         g.dicoGeneC[shName2]=GeneC(g,g1,g2)
   for g in listFeatures2:
      if g.genre=='gene':
         g.dicoGeneC[shName1]=GeneC(g,g2,g1)

   fileO=open(path+shName1+'.'+shName2+'.orth.pairs','r')

   for ligne in fileO:
      ligneList=re.split(" ",ligne)
      if ligneList[11]=='invert' or ligneList[11]=='-1':
         sens=-1
      else:
         sens=1
      gene1=listFeatures1[int(ligneList[4])-1]
      gene1c=gene1.dicoGeneC[shName2]
      gene2=listFeatures2[int(ligneList[9])-1]
      gene2c=gene2.dicoGeneC[shName1]
      if ligneList[2]!='00000':         # BDBH
         gene1c.setBDBH(gene2c,int(ligneList[10]),sens)
         gene2c.setBDBH(gene1c,int(ligneList[10]),sens)
      else:                        # 30%
         gene1c.addHomologues30(gene2c,int(ligneList[10]),sens)
         gene2c.addHomologues30(gene1c,int(ligneList[10]),sens)
   fileO.close()

######
## On complete l'objet GENEC en lui associant des BLOCK qu'on cree
# a partir du fichier '.gene.synt'
def syntheseBlocs2Genomes(genome1,genome2,path,outfileName=''):

   if genome1.name<genome2.name:
      g1=genome1
      g2=genome2
   else:
      g2=genome1
      g1=genome2
   Name1=g1.name
   Name2=g2.name
   listFeatures1=g1.features
   listFeatures2=g2.features

   ## blocs1 and blocs2
   blocs1=[]
   blocs2=[]

   bothCentro=g1.centromere() and g2.centromere()

   ### We fill them
   ## on cree les blocks !!!!!
   syntonFile=open(path+Name1+'.'+Name2+'.orth.synt','r')
   
   for i, ligne in enumerate(syntonFile):
      ligneList=re.split(" ",ligne)
      chromo1=int(ligneList[0])
      chromo2=int(ligneList[1])
      bloc1=[]
      bloc2=[]
      for j in range(len(ligneList)/5):
         bloc1.append(int(ligneList[5*j+3]))
         bloc2.append(int(ligneList[5*j+5]))
      min1=min(bloc1)
      max1=max(bloc1)
      min2=min(bloc2)
      max2=max(bloc2)

      #the sign of bloc 2
      sign2=sign(bloc1,bloc2,listFeatures1,listFeatures2,bothCentro)
      # un block contient des genes et des ancres
      newBloc1=Block('B'+"%.5d"%(i+1)+'G1',1,chromo1,(min1,max1))
      newBloc2=Block('B'+"%.5d"%(i+1)+'G2',sign2,chromo2,(min2,max2))
      for j in range(min1,max1+1):
         g1lF=listFeatures1[j-1]
         if g1lF.genre=='gene':
            g1lF.dicoGeneC[Name2].addBlocks(newBloc1)
            newBloc1.addGenes(g1lF.dicoGeneC[Name2])
         if j in bloc1 and g1lF not in newBloc1.ancres:
            newBloc1.addAncres(g1lF.dicoGeneC[Name2])
      for j in range(min2,max2+1):
         g2lF=listFeatures2[j-1]
         if g2lF.genre=='gene':
            g2lF.dicoGeneC[Name1].addBlocks(newBloc2)
            newBloc2.addGenes(g2lF.dicoGeneC[Name1])
         if j in bloc2 and g2lF not in newBloc2.ancres:
            newBloc2.addAncres(g2lF.dicoGeneC[Name1])

      newBloc1.genes.sort(lambda x,y:cmp(x.position,y.position))
      newBloc1.ancres.sort(lambda x,y:cmp(x.position,y.position))
      newBloc2.genes.sort(lambda x,y:cmp(x.position,y.position))
      newBloc2.ancres.sort(lambda x,y:cmp(x.position,y.position))
      # we look if it is included in another
      includ(newBloc1,blocs1)
      overlap(newBloc1,blocs1)
      includ(newBloc2,blocs2)
      overlap(newBloc2,blocs2)

      newBloc2.setBlockG(newBloc1)
      newBloc1.setBlockG(newBloc2)
      
      #newBloc1.setBornesBDBH()
      #newBloc2.setBornesBDBH()
      blocs1.append(newBloc1)
      blocs2.append(newBloc2)
   removeIncludedinBoth(blocs1,blocs2,listFeatures1,listFeatures2,Name1,Name2)
   blocs2=triBloc(blocs2)

   print "BLOCKS  ",Name1,"/",Name2
   return blocs1,blocs2

## Definition du signe d'un bloc
def sign(bloc1,bloc2,lf1,lf2,bothCentro):
   min1=min(bloc1)
   max1=max(bloc1)
   min2=min(bloc2)
   max2=max(bloc2)
   bloc1C=bloc1[:]
   bloc2C=bloc2[:]
   signPos=0
   signNeg=0
   if bothCentro:
      while min1 in bloc1C:
         i1=bloc1C.index(min1)
         if bloc2C[i1]==min2 and lf1[min1-1].strand==lf2[min2-1].strand:
            signPos+=1
         elif bloc2C[i1]==max2 and lf1[min1-1].strand!=lf2[max2-1].strand:
            signNeg+=1
         bloc1C.remove(min1)
         bloc2C.remove(bloc2C[i1])
      while max1 in bloc1C:
         i1=bloc1C.index(max1)
         if bloc2C[i1]==max2 and lf1[max1-1].strand==lf2[max2-1].strand:
            signPos+=1
         elif bloc2C[i1]==min2 and lf1[max1-1].strand!=lf2[min2-1].strand:
            signNeg+=1
         bloc1C.remove(max1)
         bloc2C.remove(bloc2C[i1])
   else:
      while min1 in bloc1C:
         i1=bloc1C.index(min1)
         if bloc2C[i1]==min2:
            signPos+=1
         elif bloc2C[i1]==max2:
            signNeg+=1
         bloc1C.remove(min1)
         bloc2C.remove(bloc2C[i1])
      while max1 in bloc1C:
         i1=bloc1C.index(max1)
         if bloc2C[i1]==max2:
            signPos+=1
         elif bloc2C[i1]==min2:
            signNeg+=1
         bloc1C.remove(max1)
         bloc2C.remove(bloc2C[i1])
   if signPos<signNeg:
      return -1
   elif signPos>signNeg:
      return 1
   else:
      return 0

## Mise a jour des inclusion
def includ(bl,blocs):
   for bloc in blocs:
      if (int(bloc.bornes[0])<=int(bl.bornes[0])
         and int(bloc.bornes[1])>=int(bl.bornes[1])):
         bloc.addIncluded(bl) 
         bl.addIncluding(bloc)
      elif (int(bloc.bornes[0])>=int(bl.bornes[0])
           and int(bloc.bornes[1])<=int(bl.bornes[1])):
         bl.addIncluded(bloc)  
         bloc.addIncluding(bl)
## Mise a jour des overlap
def overlap(bl,blocs):
   if bl.including==[]:
      for bloc in blocs:
         if bloc.including==[]:
            if (int(bloc.bornes[0])>int(bl.bornes[0])
               and int(bloc.bornes[0])<=int(bl.bornes[1])
               and int(bloc.bornes[1])>int(bl.bornes[1])):
               bloc.setOverlapL(bl)
               bl.setOverlapR(bloc)
            elif (int(bloc.bornes[0])<int(bl.bornes[0])
                 and int(bloc.bornes[1])<int(bl.bornes[1])
                 and int(bloc.bornes[1])>=int(bl.bornes[0])):
               bloc.setOverlapR(bl)
               bl.setOverlapL(bloc)

## Supprime le bloc si il est inclus ds les deux genomes
# partageant les memes ancres que le blocs qui l'inclue 
def removeIncludedinBoth(blocs1,blocs2,listFeatures1,listFeatures2,Name1,Name2):
   for bloc in blocs1[:]:
      if bloc.including!=[]:
         d1=0
         for bl in bloc.including:
            if (bloc.bornes!=bl.bornes
               and listInclus(bloc.ancres,bl.ancres)
               and d1==0):
               bloc2=bloc.blockG
               if bloc2.including!=[]:
                  d2=0
                  for bl2 in bloc2.including:
                     if (bloc2.bornes!=bl2.bornes
                        and listInclus(bloc2.ancres,bl2.ancres)
                        and d2==0):
                        for x in range(bloc.bornes[0],bloc.bornes[1]+1):
                           g=listFeatures1[x-1]
                           if g.genre=='gene':
                              g.dicoGeneC[Name2].removeBlocks(bloc)
                        for x in range(bloc2.bornes[0],bloc2.bornes[1]+1):
                           g=listFeatures2[x-1]
                           if g.genre=='gene':
                              g.dicoGeneC[Name1].removeBlocks(bloc2)
                        removeB(bloc)
                        removeB(bloc2)
                        blocs1.remove(bloc)
                        blocs2.remove(bloc2)
                        d1=1
                        d2=1
def listInclus(liste1,liste2): # liste1 inclus ds liste2
   for bl in liste1:
      if bl not in liste2:
         return 0
   return 1
def removeB(bloc): # le retire ds ses relations possibles avec d'autres blocs
   for bl in bloc.including:
      bl.removeIncluded(bloc)
   for bl in bloc.included:
      bl.removeIncluding(bloc)
   if bloc.overlapR!=0:
      bloc.overlapR.setOverlapL(0)
   if bloc.overlapL!=0:
      bloc.overlapL.setOverlapR(0)
## Supprime le bloc si il est inclus ds les deux genomes
# partageant les memes ancres que le blocs qui l'inclue


## Remets les blocs ds l'ordre
def triBloc(listBlocs):
   blocsTrie=[listBlocs[0]]
   listBlocs.remove(listBlocs[0])
   for bloc in listBlocs:
      p=0
      if bloc.bornes[0]<blocsTrie[0].bornes[0]:
         blocsTrie.insert(0,bloc)
      else:
         for i in range(len(blocsTrie)-1):
            if (bloc.bornes[0]>=blocsTrie[i].bornes[0]
               and bloc.bornes[0]<blocsTrie[i+1].bornes[0]):
               blocsTrie.insert(i+1,bloc)
               p=1
               break
         if bloc.bornes[0]>=blocsTrie[-1].bornes[0] and not p:
            blocsTrie.append(bloc)
   return blocsTrie

# on supprime les blocs sans aucun homo ds G3
def suppressionBlocs(blocs1,blocs2,genomes3,outfileName):
   out=open(outfileName,'a')
   out.write('Deleted blocs (because they contain genes not present in outgroups)\n')
   out.write('previous nb of blocks '+str(len(blocs1))+'\n')
   for b in blocs1[:]:
      bG=b.blockG
      bAncres=b.ancres
      bGAncres=bG.ancres
      res=1
      for gc in bAncres:
         for genome in genomes3:
            gc3=gc.gene.dicoGeneC[genome.name]
            if gc3.blocks!=[]:
               res=0
               break
         if not res:
            break
      if res:
         for gc in bGAncres:
            for genome in genomes3:
               gc3=gc.gene.dicoGeneC[genome.name]
               if gc3.blocks!=[]:
                  res=0
                  break
            if not res:
               break
      if res:
         out.write(b.__str__()+bG.__str__()+'\n')
         blocs1.remove(b)
         removeB(b)
         blocs2.remove(bG)
         removeB(bG)
   out.write('new nb of blocks '+str(len(blocs1))+'\n')
   out.close()
   return (blocs1,blocs2)

def nbDifferentPack(blocs1,blocs2,out):
   out.write('Overlapping\tUnsigned\tIncluding\tIncluded\tIncludedG\tTelomeric\tBoth\tBothG\n')
   res=[0,0,0,0,0,0,0,0]
   for b in blocs1+blocs2:
      if b.overlapL!=0 or b.overlapR!=0:
         res[0]+=1
      if b.sign==0:
         res[1]+=1
      if b.included!=[]:
         res[2]+=1
      if b.including!=[] and b.isInSubtelo():
         res[6]+=1
      elif b.isIncludedG() and b.isInSubtelo():
         res[7]+=1
      elif b.including!=[]:
         res[3]+=1
      elif b.isIncludedG():
         res[4]+=1
      elif b.isInSubtelo():
         res[5]+=1
   print res
   out.write('\t'.join(map(str,res))+'\n')



















