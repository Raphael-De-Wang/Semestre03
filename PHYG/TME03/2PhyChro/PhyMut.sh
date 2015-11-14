#!/bin/bash

## we compute gene families
# $1:CladeName $2:Delta $3:minPercentOfSimi
./InterOrtho.py $1 $2 $3

cd BioTools/
echo `ls`
datadir=../../../$1/12OrthFamilies/Delta$2

## we aligne protein sequences of each family
for f in $datadir/Prot_* ; do
   Muscle3.8.31/./muscle3.8.31_i86linux64 -in $f -out ${f/%.fasta/.muscle}
#   Muscle3.8.31/./muscle3.8.31_i86darwin64 -in $f -out ${f/%.fasta/.muscle}
done

## we clean the alignments
for f in $datadir/*.muscle ; do
   Gblocks_Linux_0.91b/./Gblocks $f $datadir/
#   Gblocks_OSX_0.91b/./Gblocks $f $datadir/
done

## we concatenate the alignments
cat  $datadir/*-gb > $datadir/AlignConcat
./concat.py $datadir/AlignConcat

## we reconstruct the tree with phyML
PhyML_3.0/./PhyML_3.0_linux64 -i $datadir/AlignConcat.fasta -b 100 -d aa -s best -q
#PhyML_3.0/./PhyML_3.0_macOS_i386 -i $datadir/AlignConcat.fasta -b 100 -d aa -s best -q

exit
