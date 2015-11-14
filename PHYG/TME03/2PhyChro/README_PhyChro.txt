

     PhyChro  -   January 2015


-------------------------------------------------------------------------------
Table of Contents
-------------------------------------------------------------------------------
 1. Overview 
 2. Input Description
 3. PhyChro Usage 
 4. Output Description
 5. Reconstruction using PhyML using Families of Syntenic Homologs 


-------------------------------------------------------------------------------
 1. Overview
-------------------------------------------------------------------------------
PhyChro is a software tool that reconstructs the phylogenetic tree associated 
to n species using the information contained in their synteny blocks 
(reconstructed for each pairwise comparisons of these n genomes).

PhyChro requires genomes description files and several sets of synteny blocks
(one per each pairwise comparison). 
(i)   If synteny blocks have been computed with SynChro, you can directly run 
      PhyChro (see 3.).
(ii)  To reconstruct synteny blocks with SynChro, look at the SynChro's webpage 
      (http://www.lgm.upmc.fr/CHROnicle/SynChro.html).
(iii) To reconstruct the tree from others synteny blocks, you will need to 
      convert them into a specific format (see 2. for the format of input 
      files). 

PhyChro is fast enough to be used, on a desk computer, on large eukaryotic 
genomes (15 minutes for 13 vertebrates and 8 minutes for 21 yeasts) but it is 
very ram-consuming and therefore the number of compared species cannot be 
too high (depending on the memory space available on the computer).

The output tree is in the newick format 
(http://en.wikipedia.org/wiki/Newick_format, see 4.)


-------------------------------------------------------------------------------
 2. Input Description
-------------------------------------------------------------------------------
Each species name is abbreviated by a 4 letter short name ("NAME" in the 
following)

Three different files are needed to run PhyChro :

(1) The 'NAME.def' files (one file for each genome)
-> it corresponds to the list of the genetic features (one per line) of the 
   genome NAME 
for each feature, a line with 10 columns (tab-separated values):
   1 type      gene (protein coding genes),tRNA,ncRNA,rRNA,repeatR,
               LTR,pseudoG,centromere 
   2 name      feature name
   3 chr       chromosome (or scaffold) number or name
               (which can start with 000 or 001)
   4 start     start position in nucleotide, along the chromosome (or 
               scaffold)
   5 end       end position in nucleotide, along the chromosome (or 
               scaffold)
               (start and end are optional, they can be equal to 0)
   6 strand    + or - 
               (useful for the definition of the block sign;
               optional, they can be equal to 0)
   7 sens      f (for "from" centromere) or t (for "toward" centro)
               (useful for the orientation (red and green color) in 
               the Synteny map; optional, they can be equal to 0)
   8 IDg/chr   protein-coding gene number in the chromosome (or 
               scaffold) 
   9 IDg/all   protein-coding gene number in the genome
   10 IDf/all  feature number in the genome

  (-- the 5 first lines of CHROnicle/Yeast/01Genomes/YALI.def:
type	name		chr	start	end	strand	sens	IDg/chr	IDg/all	IDf/all
gene	YALI0A00110g	001	2659	5277	+	t	00001	00001	00001
gene	YALI0A00132g	001	7045	8880	+	t	00002	00002	00002
gene	YALI0A00154g	001	11559	12653	+	t	00003	00003	00003
gene	YALI0A00176g	001	15861	18419	+	t	00004	00004	00004 --)


(2) The 'NAM1.NAM2.orth.pairs' files (one file for each pairwise genome 
                                                                comparison)
-> it corresponds to the pairs (one per line) of homologous genes conserved in 
   synteny, between the two genomes NAM1 and NAM2
for each pair, a line with 12 columns (space-separated values):
   1 name1     gene name of the gene of genome NAM1
   2 chr1      chromosome (or scaffold) number in the genome NAM1
   3 IDr/all1  rbh number in the genome NAM1
   4 IDg/all1  gene number in the genome NAM1
   5 IDf/all1  feature number in the genome NAM1
   6,7,8,9,10  the same for its homolog in NAM2
   11 %simi    percentage of similarity
   12 sstr     sameW or inverted (same strand or not)

  (-- the 3 first lines of CHROnicle/Yeast/11Blocks/Delta3/OrthBlocks/CAAL.CADU.orth.pairs:
CAL0005939 001 00001 00006 00006 CD36_00070 001 00003 00005 00006 92 sameW
CAL0005938 001 00002 00007 00007 CD36_00080 001 00004 00006 00007 91 sameW
CAL0005936 001 00003 00008 00008 CD36_00090 001 00005 00007 00008 87 sameW --)

(3) The 'NAM1.NAM2.orth.synt' files (one file for each pairwise genome
                                                                comparison)
-> it corresponds to the description of the synteny blocks (one per line)
for each block, a line with 2+5*#anchors columns (space-separated values):
   1 chr1       chromosome (or scaffold) number 1
   2 chr2       chromosome (or scaffold) number 2
   3 name1      gene name 1 of the first anchor pair
   4 IDf/all1   associated IDf/all
   5 name2      gene name 2 of the first anchor pair
   6 IDf/all2   associated IDf/all
   7 %simi      percentage of similarity
   8,9,10,12,13 the second pairs of anchors defining the block
   ....         other pairs if they exist

  (-- the 3 first lines of CHROnicle/Yeast/11Blocks/Delta3/OrthBlocks/CAAL.CAGL.orth.synt:
001 003 CAL0005939 00006 CAGL0C03608g 00587 64 CAL0005938 00007 CAGL0C03630g 00588 76
001 009 CAL0005191 00033 CAGL0I03388g 02749 50 CAL0005189 00034 CAGL0I03366g 02748 64
001 008 CAL0005165 00042 CAGL0H05841g 02388 49 CAL0005161 00045 CAGL0H05797g 02386 71 CAL0005160 00046 CAGL0H05885g 02390 52 --)

These files can be directly created from GenBank/EMBL/Fasta files, look at
the "CHROnicle/Programs/0Convert2inputF/README_Convert.txt" file for a
detailed description of how to generate them.


-------------------------------------------------------------------------------
 3. Usage 
------------------------------------------------------------------------------- 
If you use PhyChro using as input synteny blocks computed with SynChro:
Compute the phylogenetic tree associated to your species (those with a 
.def description file present in the "CHROnicle/CladeName/01Genomes/" 
directory), by going to the "CHROnicle/Programs/2PhyChro/" directory and 
executing:

   ./PhyChro.py 'CladeName' Delta 'outputName' a

Delta:  is a int > 0 associated to the reconstruction of your synteny blocks 
        that you wish to use for your phylogenetic reconstruction
    a:  is an integer equal to 0 or 1 
        that allows you to reconstruct the whole tree (0) 
        or only the tree associated to a subset of genomes (1)

This will create two output files named 'outputName.outtree' and 
'outputName.out' in the "CHROnicle/CladeName/20Trees/Delta/" directory.

(ex: ./PhyChro.py Yeast 3 tree1 0)

------------------------------------------------------------------------------- 
If you use PhyChro with other kinds of synteny blocks:
Compute the phylogenetic tree associated to all your species (i.e. those 
with their .def, .orth.pairs and .orth.synt description files present in a 
"Path" directory), by going to the "CHROnicle/Programs/2PhyChro/" directory and 
executing:

   ./PhyChro.py 'Path' 'outputName' a

    a: is a int equal to 0 or 1 
       that allows you to reconstruct the whole tree (0) 
       or only the tree associated to a subset of genomes (1)

This will create a two output files 'outputName.outtree' and 
'outputName.out' in the "Path" directory.

(ex: ./PhyChro.py ./Data/ tree1 0)


-------------------------------------------------------------------------------
 4. Output Description
-------------------------------------------------------------------------------
In the "CHROnicle/CladeName/20Trees/Delta" (or "Path") directory

(1) tree1.outtree
-> corresponds to the newick tree (with branch length and confidence score)

(2) tree1.out
-> summarizes the successive steps of the reconstruction
showing for each step: 
   - the finc, fcomp, (finc+1)/(fcomp+1) values associated to each  
     pairwise genome comparison
   - the list of the n/2 pairs of genomes with the smallest (finc+1)/(fcomp+1) 
     values sorted according to their finc values


-------------------------------------------------------------------------------
 5. Reconstruction using PhyML using Families of Syntenic Homologs (PhyMut)
-------------------------------------------------------------------------------
To compute the PhyML phylogenetic tree associated to all your species and based 
on families of syntenic homologs deduced from synteny blocks (reconstructed 
with Delta), go to "CHROnicle/Programs/2PhyChro/" and execute:

   ./PhyMut.sh CladeName Delta MiniPercentOfSimi

where   MiniPercentOfSimi is a number between 40 and 100, corresponding to the 
        minimum of average similarity required between orthologs of a same 
        family

(ex:    ./PhyMut.sh Yeast 3 80)

The associated result files (the syntenic homolog families, the cleaned 
alignments, the PhyML outputs) are in 
"CHROnicle/CladeName/12OrthFamilies/Delta?/"

(Notice that some programs, as Gblocks, used by PhyMut.sh are compiled to run 
on a 32 bits computer, you may need to install ia32-libs)

