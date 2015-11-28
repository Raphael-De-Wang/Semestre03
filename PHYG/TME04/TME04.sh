#!env bash

alias archaeopteryx="java -jar /users/nfs/Etu9/3404759/Workspace/toolkits/Archaeopteryx/forester_1038.jar"

align () {

    python Exercise02.py -if $1 -is TME4bis/RB.list -of RB.seq
    clustalo -i RB.seq --seqtype=Protein --outfile aln-phylip --outfmt=phy --force

}

protdist () {

    phylip protdist <<EOF
$1
Y
EOF

    mv outfile protDist_outfile

}

neighbor () {

    phylip neighbor <<EOF
$1
Y
EOF

    mv outfile NJ_outfile
    mv outtree NJ_outtree

}

Q2 () {
    align TME4bis/RB_sequences/PF01599.12.seq
    protdist aln-phylip
    neighbor protDist_outfile
    archaeopteryx NJ_outtree
}

Q3 () {

    for f in TME4bis/RB_sequences/PF*seq ; do 
	fname=(`echo $f | tr '/' ' '`)
	pfam=(`echo ${fname[2]} | tr '.' ' '`)
	align $f && mv aln-phylip ${pfam[0]}.phylip
    done

    export HTTP_PROXY="http://proxy:3128"

    url=`python species_heatmap.py -if PF*phylip -is TME4bis/RB.list -of heatmap`

    firefox $url

    url=`python species_heatmap.py -if PF*phylip -is TME4bis/RB.list -of heatmap --filter 0.9 0.8`

    firefox $url

    python species_heatmap.py -if PF*phylip -is TME4bis/RB.list -of heatmap --filter 0.9 0.8 --no-heatmap --concat aln-phylip

    protdist aln-phylip
    neighbor protDist_outfile
    archaeopteryx NJ_outtree
}

Q4 () {
    
}

# Q2 
# Q3 


