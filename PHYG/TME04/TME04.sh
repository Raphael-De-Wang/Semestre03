#!env bash

alias archaeopteryx="java -jar /users/nfs/Etu9/3404759/Workspace/toolkits/Archaeopteryx/forester_1038.jar"

align () {

    python Exercise02.py -if $1 -is TME4bis/RB.list -of RB.seq
    clustalo -i RB.seq --seqtype=Protein --outfile $2 --outfmt=phy --force
    rm -f RB.seq
}

alignMulti () {
    python Exercise04.py -if `ls $1/*` -of RB.seq
    clustalo -i RB.seq --seqtype=Protein --outfile $2 --outfmt=phy --force
    rm -f RB.seq
}

protdist () {

    phylip protdist <<EOF
$1
Y
EOF

    mv outfile $2

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
    align TME4bis/RB_sequences/PF01599.12.seq aln-phylip
    protdist aln-phylip protDist_outfile
    neighbor protDist_outfile
    archaeopteryx NJ_outtree
    rm -f aln-phylip
}

Q3 () {

    phylip_path="TME4bis/RB_sequences_phylip"
    mkdir -p $phylip_path

    for f in TME4bis/RB_sequences/PF*seq ; do 
	fname=(`echo $f | tr '/' ' '`)
	pfam=(`echo ${fname[2]} | tr '.' ' '`)
	align $f $phylip_path/${pfam[0]}.phylip
    done

    export HTTP_PROXY="http://proxy:3128"

    url=`python species_heatmap.py -if $phylip_path/PF*phylip -is TME4bis/RB.list -of heatmap`

    firefox $url

    url=`python species_heatmap.py -if $phylip_path/PF*phylip -is TME4bis/RB.list -of heatmap --filter 0.9 0.8`

    firefox $url

    python species_heatmap.py -if $phylip_path/PF*phylip -is TME4bis/RB.list -of heatmap --filter 0.9 0.8 --no-heatmap --concat aln-phylip

    protdist aln-phylip protDist_outfile
    neighbor protDist_outfile
    archaeopteryx NJ_outtree
    rm -f aln-phylip
}

Q4 () {
    pfam_path='TME4bis/RB_clades/'
    # for pfam in `ls $pfam_path` ; do 
    # alignMulti "$pfam_path$pfam" "$pfam_path$pfam/$pfam.phylip"
    # done

    pfam_list=(`ls $pfam_path`)
    pfam_pair_path='TME4bis/RB_clades_pairs/'
    pfam_pair_prodist_path='TME4bis/RB_clades_pairs_prodist/'
    pfam_pair_tree_path='TME4bis/RB_clades_pairs_tree/'

    mkdir -p $pfam_pair_path $pfam_pair_prodist_path $pfam_pair_tree_path

    for p1 in `seq 0 $(expr ${#pfam_list[@]} - 1)` ; do
	for p2 in `seq $(expr $p1 + 1) $(expr ${#pfam_list[@]} - 1)` ; do
	    ph1="$pfam_path${pfam_list[$p1]}/${pfam_list[$p1]}.phylip"
	    ph2="$pfam_path${pfam_list[$p2]}/${pfam_list[$p2]}.phylip"
	    pfam_pair_phylip="$pfam_pair_path${pfam_list[p1]}_${pfam_list[p2]}.phylip"
	    pfam_pair_prodist="$pfam_pair_prodist_path${pfam_list[p1]}_${pfam_list[p2]}.prodist"
	    pfam_pair_outfile="$pfam_pair_tree_path${pfam_list[p1]}_${pfam_list[p2]}.outfile"
	    pfam_pair_outtree="$pfam_pair_tree_path${pfam_list[p1]}_${pfam_list[p2]}.outtree"
	    # echo "python concat_phylip.py -if $ph1 $ph2 -of $pfam_pair_phylip"
	    # echo "protdist $pfam_pair_phylip $pfam_pair_prodist"
	    # echo "neighbor $pfam_pair_prodist && mv NJ_outfile $pfam_pair_tree"
	    # python concat_phylip.py -if $ph1 $ph2 -of $pfam_pair_phylip
	    # protdist $pfam_pair_phylip $pfam_pair_prodist
	    # neighbor $pfam_pair_prodist && mv NJ_outfile $pfam_pair_outfile && mv NJ_outtree $pfam_pair_outtree 
	    python concat_phylip.py -if $ph1 $ph2 -of $pfam_pair_phylip && protdist $pfam_pair_phylip $pfam_pair_prodist && neighbor $pfam_pair_prodist && mv NJ_outfile $pfam_pair_outfile && mv NJ_outtree $pfam_pair_outtree
	done
    done

}

# Q2 
# Q3 
Q4

