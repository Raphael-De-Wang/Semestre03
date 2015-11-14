#!/usr/bin/env python
# -*- coding: utf-8 -*-
## n06CalculScore
## Compare the linkedBlocs found for 2 genomes to a outgroup


def brScore(gaucheBl,rg,droitBl,rd,genome3,macro,blocs=[]):

   # blocs auquelles appartiennent les genes de G et D in G1/G3
   blsG=[]
   blsD=[]
   
   if 'v' in gaucheBl.name:
      gaucheBl=gaucheBl.block
   if 'v' in droitBl.name:
      droitBl=droitBl.block
       
   gaucheG=gaucheBl.ancres[:]
   droitG=droitBl.ancres[:]

   ## bloc droit inclus ds gauche
   if droitBl in gaucheBl.included:
      for g in gaucheG:
         if g not in droitG:
            gc=g.gene.dicoGeneC[genome3.name]
            blsG+=[x for x in gc.blocks if x not in blsG]
      for g in droitG:
         gc=g.gene.dicoGeneC[genome3.name]
         for b in gc.blocks:
            if gc in b.ancres and b not in blsD:
               blsD+=[b]
      if [x for x in blsD if x in blsG]!=[]:#intersection non nul
         return 0.4,genome3.name
      else:
         return 0,genome3.name
   ## bloc gauche inclus ds droit  
   if gaucheBl in droitBl.included:
      for g in droitG:
         if g not in gaucheG:
            gc=g.gene.dicoGeneC[genome3.name]
            blsD+=[x for x in gc.blocks if x not in blsD]
      for g in gaucheG:
         gc=g.gene.dicoGeneC[genome3.name]
         for b in gc.blocks:
            if gc in b.ancres and b not in blsG:
               blsG+=[b]
      if [x for x in blsG if x in blsD]!=[]:#intersection non nul
         return 0.4,genome3.name
      else:
         return 0,genome3.name

   ## bloc non inclus
   if not rg:
      #print 'CalculScore 1', br1.gauche.name
      gaucheG.reverse()
   if rd:
      #print 'CalculScore 2', br1.droit.name
      droitG.reverse()
   
   ## si chevauchement
   if gaucheBl.bornes[0]<droitBl.bornes[0] and gaucheBl.bornes[1]>=droitBl.bornes[0]:
      ch=(droitBl.bornes[0],gaucheBl.bornes[1])
   elif gaucheBl.bornes[0]>droitBl.bornes[0] and gaucheBl.bornes[0]<=droitBl.bornes[1]:
      ch=(gaucheBl.bornes[0],droitBl.bornes[1])
   else:
      ch=0
       
   l1=l2=r1=r2=0
   for g in gaucheG:
      if (ch!=0 and (g.position<ch[0] or g.position>ch[1])) or ch==0:
         gc=g.gene.dicoGeneC[genome3.name]
         ## on n'a pas encore atteint une 2eme ancre appartenant a un bloc ds l'outgroup
         if not l2:
            ## on n'a pas encore atteint une 1ere ancre
            if not l1:
               Blsl1,insideSDicl1=milieuORnot(gc,(1-rg))
               if Blsl1!=[]:
                  l1=gc
            else:
               Blsl2,insideSDicl2=milieuORnot(gc,(1-rg))
               if Blsl2!=[]:
                  l2=gc
   for g in droitG:
      if (ch!=0 and (g.position<ch[0] or g.position>ch[1])) or ch==0:
         gc=g.gene.dicoGeneC[genome3.name]
         ## on n'a pas encore atteint une 2eme ancre appartenant a un bloc ds l'outgroup
         if not r2:
            ## on n'a pas encore atteint une 1ere ancre
            if not r1:
               Blsr1,insideSDicr1=milieuORnot(gc,rd)
               if Blsr1!=[]:
                  r1=gc
            else:
               Blsr2,insideSDicr2=milieuORnot(gc,rd)
               if Blsr2!=[]:
                  r2=gc
   if not l1:
      Blsl1,insideSDicl1=milieuORnot(l1,(1-rg))
   if not l2:  
      Blsl2,insideSDicl2=milieuORnot(l2,(1-rg))
   if not r1:
      Blsr1,insideSDicr1=milieuORnot(r1,rd)
   if not r2:
      Blsr2,insideSDicr2=milieuORnot(r2,rd)
   Blsl=Blsl1+[x for x in Blsl2 if x not in Blsl1]
   Blsr=Blsr1+[x for x in Blsr2 if x not in Blsr1]
   interlr=[x for x in Blsl if x in Blsr]
   res=0
   for b in interlr:
      if b not in insideSDicl2.keys():
         insideSDicl2[b]=0
      if b not in insideSDicl1.keys():
         insideSDicl1[b]=0
      if b not in insideSDicr2.keys():
         insideSDicr2[b]=0
      if b not in insideSDicr1.keys():
         insideSDicr1[b]=0
      if b not in insideSDicl1.keys() or b not in insideSDicr1.keys():
         print b,l2.position,l1.position,r1.position,r2.position,l1.genome.name,l1.genomeC.name
           
      res=max(res,cScore(insideSDicl1[b],insideSDicl2[b],insideSDicr1[b],
                          insideSDicr2[b],l1.position,1-rg,r1.position,rd,b))
   if res!=0:
      return res,genome3.name
   elif blocs!=[]:
      ##intersection nulle ms cote a cote ds G1G3 
      return exist(Blsl,Blsr,blocs),genome3.name
   elif not macro:
      #print Ancres2BlsG1+Ancre1BlsG1,Ancres2BlsD1+Ancre1BlsD1
      return 1-distanceHomos(l1,1-rg,Blsl,r1,rd,Blsr),genome3.name
   else:
      return 0,genome3.name

def cScore(iSl1,iSl2,iSr1,iSr2,l1pos,rg,r1pos,rd,b):
   if iSl1+iSl2>=4 and iSr1+iSr2>=4:
      return 1-distanceHomosBloc(l1pos,rg,r1pos,rd,b)
   elif ((iSl1+iSl2>=4 and iSr1+iSr2==3) or
        (iSl1+iSl2==3 and iSr1+iSr2>=4)):
      return 0.9-distanceHomosBloc(l1pos,rg,r1pos,rd,b)
   elif iSl1+iSl2==3 and iSr1+iSr2==3:
      return 0.8-distanceHomosBloc(l1pos,rg,r1pos,rd,b)
   elif ((iSl1+iSl2>=3 and iSr1+iSr2<=2 and iSr1+iSr2>=1) or
        (iSl1+iSl2>=1 and iSl1+iSl2<=2 and iSr1+iSr2>=3)):
      return 0.7-distanceHomosBloc(l1pos,rg,r1pos,rd,b)
   elif (iSl1+iSl2>=1 and iSl1+iSl2<=2 and
        iSr1+iSr2>=1 and iSr1+iSr2<=2):
      return 0.6-distanceHomosBloc(l1pos,rg,r1pos,rd,b)
   else:
      print 'erreur',iSl1,iSl2,iSr1,iSr2


## check the existence de (g,d) in G3G1  (ds G1 forcement � cote, not in G3)
def exist(blsG,blsD,blocsList):
   blsDG=[]
   for d in blsD:
      blsDG.append(d.blockG)
   for g in blsG:
      g3=g.blockG
      i=blocsList.index(g3)
      if g3.sign*g.sign==1 and i!=len(blocsList)-1:
         nextB=blocsList[i+1]
         if g3.chromo==nextB.chromo:
            if nextB in blsDG:
               return 0.5
      elif g3.sign*g.sign==-1 and i!=0:
         prevB=blocsList[i-1]
         if g3.chromo==prevB.chromo:
            if prevB in blsDG:
               return 0.5
   return 0



## renvoie la distance entre 4 genes (between 0 and 0.09)
def distanceHomosBloc(pos1,rg,pos2,rd,bloc3):
   posGene1,posOrtho1=posGO(pos1,bloc3,rg)
   posGene2,posOrtho2=posGO(pos2,bloc3,rd)
   dist1=abs(posOrtho2[0]-posOrtho1[0])-1
   for p1 in posOrtho1:
      for p2 in posOrtho2:
         dist1=min(dist1,abs(p2-p1)-1)
   dist2=abs(posGene2-posGene1)-1
   dist=dist1+dist2
   if dist1<=0 and dist2<=0:
      return 0
   elif dist1<=1 and dist2<=1:
      return 0.01
   elif dist1<=2 and dist2<=2:
      return 0.02
   elif dist1<=3 and dist2<=3:
      return 0.03
   elif dist1<=4 and dist2<=4:
      return 0.04
   elif dist1<=5 and dist2<=5:
      return 0.05
   else:
      return 0.09
#   if dist<=8: #max 9 gene:
#      return 0.01*dist
#   else:
#      return 0.09

## renvoie la distance entre 4 genes (between 0 and 0.09)
def distanceHomos(pos1,rg,blocs1,pos2,rd,blocs2):
   g1s=[]
   for b in blocs1:
      g1s.extend(posGO(pos1,b,rg)[1])
   g2s=[]
   for b in blocs2:
      g2s.extend(posGO(pos2,b,rd)[1])
   res=[]
   for g1 in g1s:
      for g2 in g2s:
         res.append(abs(g1-g2))
   if res==[]:
      return 1
   else:
      dist=min(res)-1
      if dist<=10: #max 10 gene:
         return 0.01*dist
      else:
         return 1
def posGO(pos1,bloc3,rg):
   g1=0
   for i in range(len(bloc3.ancres)-1):
      if (rg and bloc3.ancres[i].position<=pos1
         and bloc3.ancres[i+1].position>pos1):
         g1=bloc3.ancres[i]
      elif (not rg and bloc3.ancres[i].position<pos1
           and bloc3.ancres[i+1].position>=pos1):
         g1=bloc3.ancres[i+1]
   #print 'g1,rg,bloc3',g1,rg,bloc3
   if not g1:
      if not rg:
         g1=bloc3.ancres[0]
      else:
         g1=bloc3.ancres[len(bloc3.ancres)-1]
   posOrtho1=[]
   if g1.BDBH:
      posOrtho1.append(g1.BDBH[0].gene.posGene)
   else:
      #print g1.genome.name,g1.genomeC.name,g1.position,g1.homologues30
      for g30 in g1.homologues30:
         posOrtho1.append(g30[0].gene.posGene)
   return g1.gene.posGene,posOrtho1


## petite fonction qui permet d'acceder au Block et pas Blockv
def realB(b):
   if b.gauche==0:
      bg=0
   else:
      if 'v' in b.gauche.name:
         bg=b.gauche.block
      else:
         bg=b.gauche
   if b.droit==0:
      bd=0
   else:
      if 'v' in b.droit.name:
         bd=b.droit.block
      else:
         bd=b.droit
   return (bg,bd)

## renvoie res=[[blocs dont je suis au milieu],[a l'extremit�]]
# au milieu = a au moins 2 ancres (incluant moi meme)de l'extremit�
# g=1 pour gauche [o g ... et g=0 pour droit ... g o]
def milieuORnot(g1c,g):
   if g1c:
      res={}
      bls=g1c.blocks
      for b in bls:
         pos=[]
         for a in b.ancres:
            pos.append(a.position)
         pos.sort()
         if g:
            #print b.genes[0].genome,b.genes[0].genomeC
            if g1c.position>=pos[1]:
               res[b]=2
            else:
               res[b]=1
         else:
            #print b.genes[0].genome,b.genes[0].genomeC
            if g1c.position<=pos[-2]:
               res[b]=2
            else:
               res[b]=1
      return bls,res
   else:
      return [],{}





