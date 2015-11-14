#!/usr/bin/env python
# -*- coding: utf-8 -*-
## Script to reconstruct phylogenetic tree from synteny bloc adjacencies

import re,os,sys,math,copy
import n00Objets,n01Genomes,n01SyntheseBlocs,n06CalculScore


#### Node 
class Node(object):
   "Definition of a Node"
   nb=0
   def __init__(self,name):
      self.name=name
      self.branches=[]
      self.genomes=[] #associated to the 3 branches
      self.contradict=0
   def addBranch(self,branch):
      self.branches.append(branch)
   def addGenomes(self,genomL):
      self.genomes.append(genomL)
   def getGenomes(self):
      if len(self.name)==4:
         return [self.name]
      else:
         return (self.genomes[0],self.genomes[1])
   def getGenomes2(self):
      if len(self.name)==4:
         return [self.name]
      else:
         return self.genomes[0]+self.genomes[1]
#### Branch
class Branch(object):
   "Definition of a Branch"
   nb=0
   def __init__(self,node1,node2):
      self.length=(0,0)
      self.contradict=(0,0)
      self.genomes=[] #[[name],pair of lists] or [pair of lists,pair of lists]
      self.nodes=[node1,node2]
      node1.addBranch(self)
   def addGenomes(self,genomL):
      self.genomes.append(genomL)
   def getGenomes(self,i):
      a=self.genomes[i]
      if len(a)==1:
         return a
      else:
         return a[0]+a[1]


def ke(m):
   ## ke translate '((CAAL,CADU),CATR)' in 'CAALCADUCATR'
   #  for a easier words sorting 
   if m=='':
      return ''
   elif len(m[0])>1: ##tuple
      return ke(m[0])+ke(m[1])
   elif m[0]=='(' or m[0]==',' or m[0]==')':
      return ke(m[1:])
   else:
      return m[0]+ke(m[1:])


def _main_(argv):
   if len(argv)!=5 and len(argv)!=4:
      print "\n   4 arguments are expected to use SynCho blocks as input:\n \
      -> the cladeName      (ex: 'Yeast')\n \
      -> Delta         (ex: '3')\n \
      -> the output file name   (ex: 'Tree36')\n \
      -> 0 (to built the whole tree) / 1 (to choose which subtree)\n"
      print "\n or 3 arguments are expected to use other blocks as input:\n \
      -> the path to get the genome and the bloc files (ex: './Data/')\n \
      -> output file name           (ex: 'Tree36')\n \
      -> 0 (to built the whole tree) / 1 (to choose a subtree)\n"
      exit(1)
   if len(argv)==5:
      pathDef='../../'+argv[1]+'/01Genomes/' 
      pathPairs='../../'+argv[1]+'/11Blocks/Delta'+argv[2]+'/OrthBlocks/'
      allQ=int(argv[4])
   else:
      pathDef=argv[1]
      pathPairs=argv[1]
      allQ=int(argv[3])

   ### va chercher les differents genomes
   iniShortNameList=n01Genomes.ListShortName(pathDef)
   iniShortNameList.sort()

   if allQ!=0:   ## Question/Interface
      i=0
      for x in iniShortNameList:
         i+=1
         print i,'\t',x

      questionGenome=raw_input("Which TREE do you want to reconstruct? \
... write their numbers separated by spaces \n")
      questionLigne=re.split(" ",questionGenome)
      if len(questionLigne)<4:
         print "!!! At least 4 genomes must be compared  !!!!"
         exit(1)
      else:
         shortNameList=[]
         for num in questionLigne:
            if int(num)<=len(iniShortNameList):
               shortNameList.append(iniShortNameList[int(num)-1])
            else:
               print "!!! All number has to correspond to a genome  !!!!"
               exit(1)
   else:
      shortNameList=iniShortNameList

   shortNameList.sort()
   print 'Tree of the species', shortNameList, '\n'

   ## (1) We do the genome syntheses
   genomes=[]
   nodes={}
   for shName in shortNameList:
      g=n00Objets.Genome(shName,pathDef)
      n01SyntheseBlocs.syntheseDef(g,pathDef)
      genomes.append(g)
      nodes[shName]=Node(shName)
   l=len(genomes)
   pairsOfBlocs=[]
   for i in range(len(genomes)-1):
      pairsOfBlocs.append([])
      for j in range(i+1,len(genomes)): 
         g1=genomes[i]
         g2=genomes[j]
         n01SyntheseBlocs.synthesePairs(g1,g2,pathPairs)
         ## (2) We do the pairwise bloc syntheses 
         (blocs1,blocs2)=n01SyntheseBlocs.syntheseBlocs2Genomes(g1,g2,pathPairs)
         pairsOfBlocs[i].append((blocs1,blocs2))
         
   PIG=[]
   tPIG=[]
   nbPIG=[]
   nbtPIG=[]
   nbPIA=0
   for i in range(len(genomes)-1):
      for j in range(i+1,len(genomes)):
         outgroup=genomes[:i]
         outgroup.extend(genomes[i+1:j])
         outgroup.extend(genomes[j+1:])
         ## (3) We identify PIAs
         brkptsAssocies,brkptsList1,brkptsList2=creationBrks(pairsOfBlocs[i][j-i-1])
         nbPIA+=len(brkptsAssocies)
         ## (4) We reconstruct PIGGs
         score(genomes[i].name,genomes[j].name,brkptsAssocies,brkptsList1,
              brkptsList2,outgroup,PIG,nbPIG,tPIG,nbtPIG)
   print '\ntotal nb of PIAs\t',nbPIA
   print 'total nb of PIGGs\t',sum(nbPIG)
   print 'nb of different PIGGs\t',len(PIG),'\n\nSuccessive sister genomes:'
   
   try : 
      os.mkdir('../../'+argv[1]+'/20Trees/')
   except OSError:
      pass 
   try : 
      os.mkdir('../../'+argv[1]+'/20Trees/Delta'+argv[2]+'/')
   except OSError:
      pass
   if len(argv)==5:
      fileT=open('../../'+argv[1]+'/20Trees/Delta'+argv[2]+'/'+argv[3]+'.outtree','w')
      fileO=open('../../'+argv[1]+'/20Trees/Delta'+argv[2]+'/'+argv[3]+'.out','w')
   else:
      try :
         fileT=open('../../'+argv[1]+'/20Trees/Delta'+argv[2]+'/'+argv[1]+argv[2]+'.outtree','w')
         fileO=open('../../'+argv[1]+'/20Trees/Delta'+argv[2]+'/'+argv[1]+argv[2]+'.out','w')
      except IOError :
         fileT=open(argv[1]+'/'+argv[2]+'.outtree','w')
         fileO=open(argv[1]+'/'+argv[2]+'.out','w')
   fileO.write('(G1, G2)\t\t[finc, fcomp]\t(finc+1)/(fcomp+1)\n\n')

   ## (5) We reconstruct the tree
   uPIG=copy.deepcopy(PIG)
   unbPIG=nbPIG[:]
   finalNode=minScore(uPIG,unbPIG,nodes,0,l-3,fileO)
   fileO.close()
   
   ## (6) We compute branch lengths and confidence score
   #contradictAle(PIG,nbPIG)
   lenBcontrN(finalNode,PIG,nbPIG,tPIG,nbtPIG,1)
   arbre=printN(finalNode,1)
   print arbre,'\n'

   fileT.write(arbre+';')
   fileT.close()

## CREATION OF PIAs
def creationBrks((blocs1,blocs2)):
   brkptsList1=[]
   brkptsList2=[]
   brkptsAssocies=[] # list de couple (brk1,brk2) 
   for i,block in enumerate(blocs1):
      bk1b=0
      bk2b=0
      bk1a=0
      bk2a=0
      # bloc must be not included in G1 nor in G2 and must be signed
      if (block.including==[] and block.blockG.including==[] and
         block.blockG.sign!=0):
         #####  G1  ####
         # the previous block must do not overlap or be included
         if i!=0:
            j=i-1
            while j>=0 and (blocs1[j].overlapR==block or
                        blocs1[j].including!=[]):
               j-=1
            if j!=0:  # normal previous block
               if blocs1[j].chromo==block.chromo:
                  bk1b=n00Objets.Breakpoint(blocs1[j],block) # b for before
         # the following)block must do not overlap or be included
         if i!=len(blocs1)-1:
            j=i+1
            while j<len(blocs1) and (blocs1[j].overlapL==block or
                               blocs1[j].including!=[]):
               j+=1
            if j!=len(blocs1):  # normal next block
               if blocs1[j].chromo==block.chromo:
                  bk1a=n00Objets.Breakpoint(block,blocs1[j]) # a for after

         #####  G2  ####
         i2=blocs2.index(block.blockG)
         if i2!=0:
            j=i2-1
            while (j>=0 and 
                  (blocs2[j].overlapR==block or blocs2[j].including!=[])):
               j-=1
            if j!=0:
               if blocs2[j].chromo==block.blockG.chromo:
                  if block.blockG.sign==1:
                     bk2b=n00Objets.Breakpoint(blocs2[j],block.blockG)
                  else:
                     bk2a=n00Objets.Breakpoint(blocs2[j],block.blockG)
         if i2!=len(blocs2)-1:
            j=i2+1
            while (j<len(blocs2) and 
                  (blocs2[j].overlapL==block or blocs2[j].including!=[])):
               j+=1
            if j!=len(blocs2):  
               if blocs2[j].chromo==block.blockG.chromo:
                  if block.blockG.sign==1:
                     bk2a=n00Objets.Breakpoint(block.blockG,blocs2[j])
                  else:
                     bk2b=n00Objets.Breakpoint(block.blockG,blocs2[j])
      ## If both previous(following) blocs are defined, a PIA is defined 
      if bk1b and bk2b:
         if not sameBrkpts(bk1b,bk2b):
            if (bk1b,bk2b) not in brkptsAssocies:
               brkptsAssocies.append((bk1b,bk2b))
            if bk1b not in brkptsList1:
               brkptsList1.append(bk1b)
            if bk2b not in brkptsList2:
               brkptsList2.append(bk2b)
      if bk1a and bk2a:
         if not sameBrkpts(bk1a,bk2a):
            if (bk1a,bk2a) not in brkptsAssocies:
               brkptsAssocies.append((bk1a,bk2a))  
            if bk1a not in brkptsList1:
               brkptsList1.append(bk1a)
            if bk2a not in brkptsList2:
               brkptsList2.append(bk2a)
   return brkptsAssocies,brkptsList1,brkptsList2

## sameBrkpts check if the two breakpoints are not the same (due to a Delta cut)
def sameBrkpts(b1,b2):
   return ((b1.gauche.name[:7]==b2.droit.name[:7] and
          b1.gauche.sign*b2.droit.sign==-1 and
          b1.droit.name[:7]==b2.gauche.name[:7] and
          b1.droit.sign*b2.gauche.sign==-1) or
         (b1.gauche.name[:7]==b2.gauche.name[:7] and
          b1.gauche.sign*b2.gauche.sign==1 and
          b1.droit.name[:7]==b2.droit.name[:7] and
          b1.droit.sign*b2.droit.sign==1))

## CREATION OF PIGGs
def score(name1,name2,brkptsAssocies,brkptsList1,brkptsList2,outgroup,PIG,nbPIG,tPIG,nbtPIG):
  k=0
  score1=[]
  for b in brkptsList1:
    gaucheBl,droitBl=n06CalculScore.realB(b) # abstraction from Blockv
    scList=[]
    for i in range(len(outgroup)):
      sc,n=n06CalculScore.brScore(gaucheBl,'-' in b.gauche.name,droitBl,'-' in b.droit.name,outgroup[i],1,blocs=[])
      if sc!=0:
         scList.append((outgroup[i].name,sc))
    score1.append(scList)
  score2=[]
  for b in brkptsList2:
    gaucheBl,droitBl=n06CalculScore.realB(b) # abstraction from Blockv
    scList=[]
    for i in range(len(outgroup)):
      sc,n=n06CalculScore.brScore(gaucheBl,'-' in b.gauche.name,droitBl,'-' in b.droit.name,outgroup[i],1,blocs=[])
      if sc!=0:
         scList.append((outgroup[i].name,sc))
    score2.append(scList)
  for b1,b2 in brkptsAssocies:
    res1=score1[brkptsList1.index(b1)]
    res2=score2[brkptsList2.index(b2)]
    namel1=map(lambda x : x[0],res1) ## noise
    namel2=map(lambda x : x[0],res2) ## noise
    identikName=filter(lambda x: x in namel2,namel1) ## noise
    ## if 1 genome belong to the 2 groups of a PIGG, we do not considered it
    if identikName==[]: 
      (res1,res2)=cleanBadscore(res1,res2)
      if res1!=[] or res2!=[]:
         res1.append(name1) 
         res1.sort()
         res2.append(name2)
         res2.sort()
         if res1[0]>res2[0]:
            res1,res2=res2,res1
         if len(res1)!=1 and len(res2)!=1:
         #if (('NACA' in res1 and 'SACE' in res2) or ('NACA' in res2 and 'SACE' in res1)):
            #print 'res1,res2',res1,res2
            #print 'genome1,b1,gg1,gd1',b1.gauche.genes[0].genome.name,b1,max(map(lambda x:x.position,b1.gauche.genes)),min(map(lambda x:x.position,b1.droit.genes))
            #print 'genome2,b2,gg2,gd2',b2.gauche.genes[0].genome.name,b2,max(map(lambda x:x.position,b2.gauche.genes)),min(map(lambda x:x.position,b2.droit.genes))
            if (res1,res2) in PIG:
               nbPIG[PIG.index((res1,res2))]+=1
            else:
               PIG.append((res1,res2))
               nbPIG.append(1)
         else:
            if (res1,res2) in tPIG:
               nbtPIG[tPIG.index((res1,res2))]+=1
            else:
               tPIG.append((res1,res2))
               nbtPIG.append(1)
                

def cleanBadscore(list1,list2):
   scoreLimit=0.96
   res1=[]
   res2=[]
   for (n1,s1) in list1:
      if s1>=scoreLimit:
         res1.append(n1)
   for (n2,s2) in list2:
      if s2>=scoreLimit:
         res2.append(n2)
   return (res1,res2)

##CREATION OF THE TREE
def minScore(PIG,nbPIG,nodes,finc,l,fileO):
   ## recursive function that identify two sister genomes and update the PIGGs
   if l==0:
      nodesL=nodes.keys()
      print '\nSum of all finc =',finc
      return finish(nodes[nodesL[0]],nodes[nodesL[1]],nodes[nodesL[2]])
   else:
      ((a,b),nodes,newfinc)=bestCouple(PIG,nbPIG,nodes,fileO)
      (newPIG,newnbPIG)=PIGGsUpdate((a,b),PIG,nbPIG)
      return minScore(newPIG,newnbPIG,nodes,finc+newfinc,l-1,fileO)

## final node !
def finish(n1,n2,n3):
   n=Node(n1.name+n2.name+n3.name)
   b1=Branch(n1,n)
   b1.addGenomes(n1.getGenomes())
   b1.addGenomes((n2.getGenomes2(),n3.getGenomes2()))
   n.addGenomes(n1.getGenomes2())
   n.addBranch(b1)
   b2=Branch(n2,n)
   b2.addGenomes(n2.getGenomes())
   b2.addGenomes((n1.getGenomes2(),n3.getGenomes2()))
   n.addGenomes(n2.getGenomes2())
   n.addBranch(b2)
   b3=Branch(n3,n)
   b3.addGenomes(n3.getGenomes())
   b3.addGenomes((n1.getGenomes2(),n2.getGenomes2()))
   n.addGenomes(n3.getGenomes2())
   n.addBranch(b3)
   completion(n1)
   completion(n2)
   completion(n3)
   return n

## completion of the tree (of the arguments of nodes and branches)
def completion(n):
   if len(n.name)!=4:
      b1=n.branches[0]
      n1=b1.nodes[0]
      b2=n.branches[1]
      n2=b2.nodes[0]
      b3=n.branches[2]
      n.addGenomes(b3.getGenomes(1))
      b1.addGenomes((n2.getGenomes2(),b3.getGenomes(1)))
      b2.addGenomes((n1.getGenomes2(),b3.getGenomes(1)))
      completion(n1)
      completion(n2)

## BEST SISTER GENOMES
def bestCouple(PIG,nbPIG,nodes,fileO):
   nameList=nodes.keys()
   nbg=len(nameList)
   nameList.sort(key=ke)

   dico={}  ## {('a','b'):[dist,prox],...}
   dij={}  ## {('a','b'):nb,...}

   for i in range(nbg-1):
      for j in range(i+1,nbg):
         dico[(nameList[i],nameList[j])]=[0,0]

   for i in range(len(PIG)):
      list1=PIG[i][0]
      l1=len(list1)
      list2=PIG[i][1]
      l2=len(list2)
      for x1 in range(l1-1):
         for y1 in range(x1+1,l1):
            dico[(list1[x1],list1[y1])][1]+=nbPIG[i]
      for x2 in range(l2-1):
         for y2 in range(x2+1,l2):
            dico[(list2[x2],list2[y2])][1]+=nbPIG[i]
      for x1 in range(l1):
         for x2 in range(l2):
            if ke(list1[x1])<ke(list2[x2]):
               dico[(list1[x1],list2[x2])][0]+=nbPIG[i]
            else:
               dico[(list2[x2],list1[x1])][0]+=nbPIG[i]
   printL=dico.keys()[:]
   printL.sort(key=ke)
   for k in printL:
      dij[k]=(dico[k][0]+1.)/(dico[k][1]+1.)
      fileO.write(str(k)+'\t'+str(dico[k])+"\t%0.3f\n"%dij[k])

   printL.sort(key=lambda k: dij[k])
   smallestD=printL[:nbg/2]
   smallestD.sort(key=lambda k: dico[k][0])
   fileO.write('List of the n/2 smallest d (sorted in function of their finc value)\n'+str(smallestD)+'\n\n')

   (a,b)=smallestD[0]
   print (a,b)
   n1=nodes[a]
   n2=nodes[b]
   ab='('+a+','+b+')'
   n=Node(ab)
   n.addGenomes(n1.getGenomes2())
   n.addGenomes(n2.getGenomes2())
   b1=Branch(n1,n)
   b1.addGenomes(n1.getGenomes())
   n.addBranch(b1)
   b2=Branch(n2,n)
   b2.addGenomes(n2.getGenomes())
   n.addBranch(b2)
   del nodes[a]
   del nodes[b]
   nodes[ab]=n
   return ((a,b),nodes,dico[(a,b)][0])

## PIGGs UPDATE
def PIGGsUpdate((a,b),PIG,nbPIG):
   for i in range(len(PIG)-1,-1,-1):
      list1=PIG[i][0]
      list2=PIG[i][1]
      a1=a in list1
      a2=a in list2
      b1=b in list1
      b2=b in list2
      if a1:
         if b1:
            if len(list1)!=2:
               res=PIG.pop(i)
               res[0].remove(a)
               res[0].remove(b)
               res[0].append('('+a+','+b+')')
               res[0].sort(key=ke)
               PIG.insert(i,res)
            else:
               PIG.pop(i)
               nbPIG.pop(i)
         elif b2:
            PIG.pop(i)
            nbPIG.pop(i)
         else:
            res=PIG.pop(i)
            res[0].remove(a)
            res[0].append('('+a+','+b+')')
            res[0].sort(key=ke)
            PIG.insert(i,res)
      elif a2:
         if b2:
            if len(list2)!=2:
               res=PIG.pop(i)
               res[1].remove(a)
               res[1].remove(b)
               res[1].append('('+a+','+b+')')
               res[1].sort(key=ke)
               PIG.insert(i,res)
            else:
               PIG.pop(i)
               nbPIG.pop(i)
         elif b1:
            PIG.pop(i)
            nbPIG.pop(i)
         else:
            res=PIG.pop(i)
            res[1].remove(a)
            res[1].append('('+a+','+b+')')
            res[1].sort(key=ke)
            PIG.insert(i,res)
      elif b1:
         res=PIG.pop(i)
         res[0].remove(b)
         res[0].append('('+a+','+b+')')
         res[0].sort(key=ke)
         PIG.insert(i,res)
      elif b2:
         res=PIG.pop(i)
         res[1].remove(b)
         res[1].append('('+a+','+b+')')
         res[1].sort(key=ke)
         PIG.insert(i,res)
   return (PIG,nbPIG)

#from one node lenBs complete all branche length
def lenBcontrN(node,PIG,nbPIG,tPIG,nbtPIG,i):
   if i:
      bL=node.branches
   else:
      bL=node.branches[:2]
   for b in bL:
      if len(b.genomes[0])==2:
         b.contradict=contradictB(b,PIG,nbPIG)
         b.length=lenB(b.genomes,PIG,nbPIG)
         lenBcontrN(b.nodes[0],PIG,nbPIG,tPIG,nbtPIG,0)
      else:
         b.length=lentB(b.genomes,tPIG,nbtPIG,)

## f translate a genome NAME g into its PIGG number
def f(g,l0,l1,l2):
   if g in l0:
      return 0
   elif g in l1:
      return 1
   elif g in l2:
      return 2
   else:
      return 3

## CONFIRMATION
# for internal branches
def lenB(g,PIG,nbPIG):
   #print 'internal Branch', g[0][0],g[0][1],g[1][0],g[1][1]
   nbPIAs=0
   nbPIGGs=0
   for i,p in enumerate(PIG):
      ab=map(lambda x:f(x,g[0][0],g[0][1],g[1][0]),p[0])
      ac=map(lambda x:f(x,g[0][0],g[0][1],g[1][0]),p[1])
      if (((0 in ab and 1 in ab and 2 not in ab and 3 not in ab) and
          (0 not in ac and 1 not in ac and 2 in ac and 3 in ac)) or
         ((0 not in ab and 1 not in ab and 2 in ab and 3 in ab) and
          (0 in ac and 1 in ac and 2 not in ac and 3 not in ac))):
         nbPIAs+=nbPIG[i]*1./(len(ab)*len(ac))
         nbPIGGs+=nbPIG[i]
   #print 'nb of PIGs that confirm it', nbPIGGs
   return (nbPIAs,nbPIGGs)
# for terminal branches
def lentB(g,tPIG,nbtPIG):
#   print 'Terminal Branch', g[0]
   nbPIAs=0
   nbPIGGs=0
   for i,p in enumerate(tPIG):
      ab=map(lambda x:f(x,g[0],[],g[1][0]),p[0])
      ac=map(lambda x:f(x,g[0],[],g[1][0]),p[1])
      if (((0 in ab and 2 not in ab and 3 not in ab) and
          (0 not in ac and 2 in ac and 3 in ac)) or
         ((0 not in ab and 2 in ab and 3 in ab) and
          (0 in ac and 2 not in ac and 3 not in ac))):
         nbPIAs+=nbtPIG[i]*1./(len(ab)*len(ac))
         nbPIGGs+=nbtPIG[i]
#   print 'nb of PIGs that confirm it', nbPIGGs
   return (nbPIAs,nbPIGGs)

## CONTRADICTION
def contradictB(branch,PIG,nbPIG):
   #print 'internal Branch ', g1,g2,g3,g4
   (g1,g2)=branch.genomes[0]
   (g3,g4)=branch.genomes[1]
   nbPIGGs=0
   for i,p in enumerate(PIG):
      ab=map(lambda x:f(x,g1,g2,g3),p[0])
      ac=map(lambda x:f(x,g1,g2,g3),p[1])
      if (((2 in ab or 3 in ab) and (0 in ab or 1 in ab) and (3 in ac or 2 in ac) and (0 in ac or 1 in ac)) or
         ((0 in ab or 1 in ab) and (2 in ab or 3 in ab) and (1 in ac or 0 in ac) and (2 in ac or 3 in ac))):
         #print p[0],p[1],1./(len(ab)*len(ac)),nbPIG[i]
         nbPIGGs+=nbPIG[i]
   return nbPIGGs

#def contradictAle(PIG,nbPIG):
   #bList=[]
   #gg1,gg2,gg3,gg4=['DANR','GALG','MOND','ORYL','TAEG','TETN'],['HOMS','MACM','PANT'],['CANF','EQUC'],['MUSM','RATN']
   #bList.append(('A34',gg1,gg2,gg3,gg4))
   #bList.append(('B34',['MUSM'],['RATN'],gg3,gg1+gg2))
   #bList.append(('C34',['CANF'],['EQUC'],gg4,gg1+gg2))
   #bList.append(('D34',['HOMS','PANT'],['MACM'],gg1,gg3+gg4))
   #bList.append(('E34',['MOND'],['DANR','GALG','ORYL','TAEG','TETN'],gg2,gg3+gg4))
   #bList.append(('A24',gg1,gg3,gg2,gg4))
   #bList.append(('B24',['MUSM'],['RATN'],gg2,gg1+gg3))
   #bList.append(('C24',['CANF'],['EQUC'],gg1,gg2+gg4))
   #bList.append(('D24',['HOMS','PANT'],['MACM'],gg4,gg1+gg3))
   #bList.append(('E24',['MOND'],['DANR','GALG','ORYL','TAEG','TETN'],gg3,gg2+gg4))
   #bList.append(('A23',gg1,gg4,gg2,gg3))
   #bList.append(('B23',['MUSM'],['RATN'],gg1,gg2+gg3))
   #bList.append(('C23',['CANF'],['EQUC'],gg2,gg1+gg4))
   #bList.append(('D23',['HOMS','PANT'],['MACM'],gg3,gg1+gg4))
   #bList.append(('E23',['MOND'],['DANR','GALG','ORYL','TAEG','TETN'],gg4,gg2+gg3))
   
   #for (l,g1,g2,g3,g4) in bList:
     #nbPIGGs=0
     #for i,p in enumerate(PIG):
     #ab=map(lambda x:f(x,g1,g2,g3),p[0])
     #ac=map(lambda x:f(x,g1,g2,g3),p[1])
     #if (((2 in ab or 3 in ab) and (0 in ab or 1 in ab) and (3 in ac or 2 in ac) and (0 in ac or 1 in ac)) or
        #((0 in ab or 1 in ab) and (2 in ab or 3 in ab) and (1 in ac or 0 in ac) and (2 in ac or 3 in ac))):
        ##print p[0],p[1],1./(len(ab)*len(ac)),nbPIG[i]
        #nbPIGGs+=nbPIG[i]
     #print l,nbPIGGs
   ##print 'nb of PIGs that contradict it', contradictA
   #return nbPIGGs

## PRINT the tree in the newick format
def printN(node,i):
   if i:
      return '('+printN(node.branches[0].nodes[0],0)+','+printN(node.branches[1].nodes[0],0)+','+printN(node.branches[2].nodes[0],0)+')'
   else:
      if len(node.name)==4:
         return node.name+':'+str(node.branches[0].length[0])
      else:
         Cpig=node.branches[2].contradict
         (Ladj,Lpig)=node.branches[2].length
         return '('+printN(node.branches[0].nodes[0],0)+','+printN(node.branches[1].nodes[0],0)+')'+str((Lpig+1.)/(Lpig+Cpig+1.))+':'+str(Ladj)




_main_(sys.argv)
##./PhyChro.py ../../Yeast/01Genomes/ ../../Yeast/11Blocks/Delta3/OrthBlocks/ ../../Yeast/20Trees/Delta3/Final 0
##./PhyChro.py ../../Vertebrate/01Genomes/ ../../Vertebrate/11Blocks/Delta3/OrthBlocks/ ../../Vertebrate/20Trees/Delta3/Final 0
