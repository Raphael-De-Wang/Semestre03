#!/usr/bin/env python
# -*- coding: utf-8 -*-
## n00Objets
## Definition of Genomes, Genes, Blocks and Packs

import os,math
import n01Genomes


## Check if B is not equal to one Si
def isIncluded(bl,l):
   bornes=bl.bornes
   for i in l:
      if bornes == i.bornes:
         return i
   return 0

def intersection(liste1,liste2):
   res=[]
   for x in liste1:
      if x in liste2:
         res.append(x)
   return res


#### GENOME ####
class Genome(object):
   "Definition d'un Genome"
   nb=0
   def __init__(self,name,path):
      self.id=Block.nb  
      Block.nb+=1  
      self.name=name           ## 'KLWA'
      self.group=path          ## 'Yeast'
      (self.listChromoName,self.listChromoCentro,self.listChromoEnd)=n01Genomes.chromosomesLists(path,name)
      self.features=[] #list of "Gene"
   ## 1 if we know the centromeres position, 0 else
   def centromere(self):
      return (self.listChromoEnd!=self.listChromoCentro)
   def setFeatures(self,f):
      self.features=f
   def __str__(self):
      return self.name
   def __hash__(self):
      return self.id 

   
#### GENE ####  in fact FEATURE
class Gene(object):
   "Definition d'un Gene (ou pas)"
   nb=0
   def __init__(self,genre,name,genome,chromo,startEnd,strand,sens,position,posGene):
      self.id=Block.nb  
      Block.nb+=1
      self.genre=genre
      self.name=name            ## le "String" 'KLWA'
      self.genome=genome        ## l'objet "Genome"
      self.chromosome=chromo    ## "int" >=1 sauf si sexuel 0
      self.startEnd=startEnd    ## (1654,1867)
      self.strand=strand        ## 1 ou -1 for + or -
      self.sens=sens            ## 1 or -1 for t or f
      self.position=position    ## ds le genome
      self.posGene=posGene      ## / aux genes (important pour le score)
      ## dico des genes equivalents dans les differentes comparaisons
      self.dicoGeneC={}         ## {'KLTH':g1, "String":l'objet "GeneC",...

   def addDicoGeneC(self,genome,gc):
      self.dicoGeneC[genome]=gc
   def isSubtelo(self):
      x=self.position
      chromoList=self.genome.listChromoEnd
##      print x, self.genome.name, chromoList
      if x<=35:
         return True
      for i in chromoList:
         if (i-35<x and x<=i) or (i<x and x<=i+35):
            return True
      return False
   def distanceEnds(self):
      x=self.position
      chromoList=self.genome.listChromoEnd
      if x<=chromoList[0]:
         return min(x,chromoList[0]-x)
      for i in range(len(chromoList)-1):
         if chromoList[i]<x and x<=chromoList[i+1]:
            return min(x-chromoList[i]-1,chromoList[i+1]-x)
   def posCentro(self):  # -1 if at right of centromere or 1 (or 0)
      cE=self.genome.listChromoEnd
      cC=self.genome.listChromoCentro
      pos=self.position
      for i in range(len(cE)):
         if pos<=cE[i]:
            if pos<cC[i]:
               return -1
            elif pos==cC[i]:
               return 0
            else:
               return 1

   def __eq__(self,other):
      return self.id==other.id
   def __str__(self):
      return str(self.position)
   def __hash__(self):
      return self.id 

#### GENEC ####
class GeneC(object):
   "Definition d'un Gene dans une Comparaison 2 à 2 avec une autre genome "
   nb=0
   def __init__(self,gene,g1,g2):
      self.id=Block.nb  
      Block.nb+=1
      self.genome=g1
      self.genomeC=g2
      self.gene=gene
      self.name=gene.name
      self.startEnd=gene.startEnd
      self.position=gene.position
      self.BDBH=0             ## un gene au mieux BDBH: ("GeneC",int%,-1 or 1 for inverted or not)
      self.homologues30=[]    ## liste des homologues 30%: ("GeneC",int%,-1 or 1 for inverted or not)
##      self.homologues3G=[]  ## liste des homologues retrouve par transition (gene,-1or1)
      self.blocks=[]          ## liste

   def setBDBH(self,g,percentage,sens):
      self.BDBH=(g,percentage,sens)
   def addHomologues30(self,g,percentage,sens):
      self.homologues30.append((g,percentage,sens))
      self.homologues30.sort(key=(lambda x : x[1]),reverse=True)
   def addBlocks(self,b):
      self.blocks.append(b)
   def removeBlocks(self,b):
      self.blocks.remove(b)
   def posCentro(self):  # -1 if at right of centromere or 1 (or 0)
      return self.gene.posCentro()
   
   def __hash__(self):
      return self.id 


#### BLOCK ####
class Block(object):
   "Definition d'un Block"
   nb=0
   def __init__(self,name,sign,chromo,bornes):
      self.id=Block.nb  
      Block.nb+=1
      self.name=name        ## B0012G1 ou B0012G2
      self.sign=sign        ## soit 1 -1 0
      self.chromo=chromo
      self.bornes=bornes
      self.bornesBDBH=bornes
      self.ancres=[]        ## genes orthologues ("GeneC")
      self.genes=[]         ## ts les genes du bloc ("GeneC")
      self.included=[]      ## the blocks that are included in this block, les blocs inclus
      self.including=[]     ## the blocks including this block , les blocs l'incluant
      self.overlapR=0       ## 'Block' the number of the right (next) block overlaping this block
      self.overlapL=0       ## the number of the left (previous) block overlaping this block
      self.previous=0
      self.next=0           ## pack or block
      self.pack=0           ## "Pack"
      self.blockG=0         ## its equivalent block in the other genome
      self.blocksV=[]
      self.brkptLeft=[]     ## en supposant le block positif et si 0 son virtuel positif
      self.brkptRight=[]
      self.cyclesL=[]       ## list de left cycles: breakpoint (_,block) 
      self.cyclesR=[]
#      self.ancetre=0       ## liste des genes ds l'ordre

   def addIncluded(self,block):
      self.included.append(block)
   def removeIncluded(self,block):
      self.included.remove(block)
   def addIncluding(self,block):
      if self.including==[]:
         self.including=[block]
      elif len(self.including[0].genes)<len(block.genes):
         self.including[0].removeIncluded(self)
         self.including=[block]
      else:
         block.removeIncluded(self)
   def removeIncluding(self,block):
      self.including.remove(block)  
   def addBlocksV(self,block):
      self.blocksV.append(block)
   def addAncres(self,g):
      self.ancres.append(g)
   def addGenes(self,g):
      self.genes.append(g)
   def addBrkptLeft(self,bp):   #(_,Bl)
      self.brkptLeft.append(bp)
   def addBrkptRight(self,bp):
      self.brkptRight.append(bp)

   def setOverlapR(self,block):
      if block!=0:
         if self.overlapR!=0:
            if block.bornes[0]<self.overlapR.bornes[0]:
               self.overlapR=block
         else:
            self.overlapR=block
      else:
         self.overlapR=0
   def removeOverlapR(self):
      self.overlapR=0
   def setOverlapL(self,block):
      if block!=0:
         if self.overlapL!=0:
            if block.bornes[1]>self.overlapL.bornes[1]:
               self.overlapL=block
         else:
            self.overlapL=block
      else:
         self.overlapL=0
   def removeOverlapL(self):
      self.overlapL=0

   def setBlockG(self,bl):
      self.blockG=bl
   def setPack(self,pack):
      self.pack=pack
   def isIncludedG(self):
      res=0
      if self.blockG.included!=[]:
         res=isIncluded(self.blockG,self.blockG.included)
      return (self.blockG.including!=[] or res)  ## 1 if included 0 otherwise

   def isSubtelo(self):
      res1=False
      for g in self.genes:
         if g.gene.isSubtelo():
            res1=True
            break
      return res1
   def isInSubtelo(self):
      res1=True
      for g in self.genes:
         if not g.gene.isSubtelo():
            res1=False
            break
      res2=True
      for g in self.blockG.genes:
         if not g.gene.isSubtelo():
            res2=False
            break
      return (res1 and res2)
   
   def setPrevious(self,bl):
      self.previous=bl
   def setNext(self,bl):
      self.next=bl

   def infPos(self,other):  ## <
      if self.bornes[0]<other.bornes[0]:
         return 1
      elif self.bornes[0]>other.bornes[0]:
         return -1
      else:
         if self.bornes[1]<other.bornes[1]:
            return 1
         elif self.bornes[1]>other.bornes[1]:
            return -1
         else:
            return 0
   def centroIncluded(self): ## 1 if the block includ the centromere
      for g in self.genes:
         if g.posCentro()==0:
            return 1
      return 0

   def setBornesBDBH(self):
      j=0
      i=1
      while i:
         i=self.ancres[j]
         if i.BDBH:
            gc2=i.BDBH[0]
            if self.blockG in gc2.blocks:
               s=i.gene.position
               i=0
         j+=1
      j=len(self.ancres)-1
      i=1
      while i:  
         i=self.ancres[j]
         if i.BDBH:
            gc2=i.BDBH[0]
            if self.blockG in gc2.blocks:
               e=i.gene.position
               i=0
         j-=1
      self.bornesBDBH=(s,e)
   def __str__(self):
      res=self.name+' '+str(self.sign)+' '+str(self.chromo)+' '+str(self.bornes)
      if self.overlapR!=0:
         res=res+' OR'+self.overlapR.name
      if self.overlapL!=0:
         res=res+' OL'+self.overlapL.name
      if self.included!=[] or self.including!=[]:
         res=res+' ['
         for i in self.included:
            res=res+'['+i.name+'],'
         for i in self.including:
            res=res+i.name+','
         res=res+']'
      if self.previous!=0:
         res=res+' prev '+self.previous.name
      else:
         res=res+' prev 0'
      if self.next:
         res=res+' next '+self.next.name
      else:
         res=res+' next 0'
      return res

   def __hash__(self):
      return self.id


#### BLOCKV ####
class Blockv(object):
   "Definition d'un Block virtual" 
   def __init__(self,bl,s,o=0):  # pour le sign et o pour l'orientation des ancres
      self.id=Block.nb       ## forme 12
      Block.nb+=1
      if o==0:
         if s==0:
            self.sign=bl.sign        ## soit 1 -1
         else:
            self.sign=s
         if self.sign==1:
            self.name=bl.name+'vp'
         ## forme B0012G1v (v for virtual p or n for positif or neg)
         else:
            self.name=bl.name+'vn'

      else:
         self.sign=s
         if o==1 and s==1: 
            self.name=bl.name+'vp+'
            ## forme B0012G1v+ (+ for same orientation)
         elif o==1 and s==-1:
            self.name=bl.name+'vn+'
         elif o==-1 and s==1: 
            self.name=bl.name+'vp-'    ## forme B0012G1v- (- for inversion of the anchors order)
         elif o==-1 and s==-1:
            self.name=bl.name+'vn-'
      self.chromo=bl.chromo
      self.block=bl
      self.blockG=bl.blockG
      self.bornes=bl.bornes
      self.previous=[]
      self.next=[]

   def setSign(self,i):
      self.sign=i
   def __str__(self):  
      res=self.name+' '+str(self.sign)+' '+str(self.chromo)+' '+str(self.block.bornes)
      if self.previous!=0:
         res=res+' prev '+self.previous.name
      else:
         res=res+' prev 0'
      if self.next:
         res=res+' next '+self.next.name
      else:
         res=res+' next 0'
      return res
   def __eq__(self,other):
      if (type(other))==int:
         return 0
      else:
         return self.id==other.id
   def __hash__(self):
      return self.id


#### PACK ####
class Pack(object):
   "Definition d'un Pack" 
   def __init__(self,name,chromo,combi):
      self.id=Block.nb
      Block.nb+=1
      self.name=name        ## forme P0012G1 ou P0012G2
      self.chromo=chromo
      self.previous=0
      self.next=0
      #self.existence=existence ## 1 if telomeric, 2 if included in G2, 3 both, else 0
      self.blocks=[]   ## non virtual
      self.combinations=combi  ## list of different orders of virtual blocks
      self.scores=[] # la liste des scores respect to each combi
##        # une liste de dico (un par combi) avec le nombre de cycles par longueur(nbPath,nbCycle)
##        self.cycles=[]
      # une liste des cycles pour chaque combi
      self.cyclesCombi=[]


   def setBlocks(self,block):
      self.blocks.append(block)
      block.setPack(self)
      for bl in block.included:
         if bl not in self.blocks:
            self.blocks.append(bl)
            bl.setPack(self)
      ch=block.overlapR 
      while ch != 0:
         while ch.including!=[]:
            ch=ch.including[0]
         self.blocks.append(ch)
         ch.setPack(self)
         for bl in ch.included:
            if bl not in self.blocks:
               self.blocks.append(bl)
               bl.setPack(self)
         ch=ch.overlapR

   def insideBrkpts(self):
      listList=self.combinations
      res=[]
      for i,l in enumerate(listList):
         for j in range(len(l)-1):
            couple=(l[j],l[j+1])
            if couple in res:
               for bp in res:
                  if bp.__eq__(couple):
                     bp.addCombination(self,[i])
            else:
               newBrkpt=Breakpoint(l[j],l[j+1],self)
               newBrkpt.addCombination(self,[i])
               res.append(newBrkpt)
      return res

   def blocksLeft(self):   # on a (B,P) ou (0,P)
      listList=self.combinations
      res={}        # dico block suivit des combi!!
      for i,l in enumerate(listList):
         if l==[]:               ## on passe le block et on arrive sur le next
            if self.next==0:
               res[0]={self:[i]}
            elif self.next.name[0]=='B':
               res[self.next]={self:[i]}
            else:
               for k,v in self.next.blocksLeft().iteritems():
                  res[k]=v
                  res[k][self]=[i]
         elif l[0] in res:
            res[l[0]][self].append(i)
         else:
            res[l[0]]={}
            res[l[0]][self]=[i]
      return res

   def blocksRight(self):    # on a (P,0) ou (P,B)
      listList=self.combinations
      res={}
      for i,l in enumerate(listList):           ## on donne tout les dernier blocs ms on ne remonte pas
         if l!=[]:
            if l[-1] not in res:
               res[l[-1]]=[i]
            else:
               res[l[-1]].append(i)
      return res

   def createCycles(self):
      for i in range(len(self.combinations)):
         self.cyclesCombi.append([])


   def setScoreList(self):
      self.scores=[]
      for cl in self.cyclesCombi:
         somme=0
         som2=0
         for c in cl:
            somme+=c.length
            som2+=c.length**2
         self.scores.append(somme/len(cl)+math.sqrt((som2/len(cl))-
                                      (somme/len(cl))**2))
      #print self.name,self.scores
      
   def __str__(self):
      res=self.name+' '+str(self.chromo)
      if self.previous!=0:
         res=res+' prev '+self.previous.name
      else:
         res=res+' prev 0'
      if self.next:
         res=res+' next '+self.next.name
      else:
         res=res+' next 0'
      res=res+' bl ['
      for i in self.blocks:
         res=res+i.name+','
      res=res+']'
      res=res+' combi ['
      for i in self.combinations:
         res=res+'['
         for j in i:
            res=res+j.name+','
         res=res+']'
      return res

   def __hash__(self):
      return self.id


#### BREAKPOINT ####
class Breakpoint(object):
   "Definition d'un Breakpoint"
   
   def __init__(self,gauche,droit,pack=0):
      self.id=Block.nb
      Block.nb+=1
      self.gauche=gauche
      self.droit=droit
      self.packCombi={}
      #un dico de dicos={"Pack":{combi(entier>0):["Cycle"s]}}
      self.packCombiCycles={} 
      if pack!=0:
         self.packCombi[pack]=[]

   def addPack(self,pack):
      self.packCombi[pack]=[]
   def addCombination(self,pack,combination):
      self.packCombi[pack].extend(combination)

   def __eq__(self,other):
      if type(other)==tuple:
         if ((((self.droit!=0 and 'G1' in self.droit.name) or
            (self.gauche!=0 and 'G1' in self.gauche.name)) and
            ((other[0]!=0 and 'G2' in other[0].name) or
            (other[1]!=0 and 'G2' in other[1].name))) or
            (((self.droit!=0 and 'G2' in self.droit.name) or
            (self.gauche!=0 and 'G2' in self.gauche.name)) and
            ((other[0]!=0 and 'G1' in other[0].name) or
            (other[1]!=0 and 'G1' in other[1].name)))):
               return 0
         else:
            if self.gauche==0:
               if self.droit==0:
                  return other==(0,0)
               else:
                  if other[1]!=0:
                     return (self.droit.name==other[1].name and other[0]==0)
                  else:
                     if other!=(0,0):
                        return ((self.droit.name[:-2]==other[0].name[:-2] and self.droit.name[:-1]!=other[0].name[:-1]) and other[1]==0)
            elif self.droit==0:
               if other[0]!=0:
                  return (self.gauche.name==other[0].name and other[1]==0)
               else:
                  if other!=(0,0):
                     return ((self.gauche.name[:-2]==other[1].name[:-2] and self.gauche.name[:-1]!=other[1].name[:-1]) and other[0]==0)
            else:
               if other[0]==0 or other[1]==0:
                  return 0
               elif self.gauche.name==other[0].name and self.droit.name==other[1].name:
                  return 1
               elif self.gauche.name[-2]=='v' and self.droit.name[-2]=='v':
                  if (self.gauche.name[:-2]==other[1].name[:-2] and self.gauche.name[-1]!=other[1].name[-1]) and (self.droit.name[:-2]==other[0].name[:-2] and self.droit.name[-1]!=other[0].name[-1]):
                     return 1
      return 0
   def __str__(self):
      res='('
      if self.gauche!=0:
         res+=self.gauche.name
      else:
         res+='0'
      res+=','
      if self.droit!=0:
         res+=self.droit.name
      else:
         res+='0'
      res+=')'
      for k,v in self.packCombi.iteritems():
         res+= k.name
         for vi in v:
            res+=' '+str(vi)
         res+=' '
      return res
   
   def __hash__(self):
      return self.id




