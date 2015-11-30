#!env bash

export LD_LIBRARY_PATH=~/GTL/lib/:LD_LIBRARY_PATH

alias archaeopteryx="java -jar /users/nfs/Etu9/3404759/Workspace/toolkits/Archaeopteryx/forester_1038.jar"
alias supertree="/users/nfs/Etu9/3404759/Workspace/Semestre03/PHYG/TME04/TME4bis/treePack/superTree/supertree"

align () {

    python Exercise02.py -if $1 -is TME4bis/RB.list -of RB.seq
    clustalo -i RB.seq --seqtype=Protein --outfile $2 --outfmt=phy --force
    rm -f RB.seq
}

align_dir () {
    python Exercise04.py -if `ls $1/*` -of RB.seq
    clustalo -i RB.seq --seqtype=Protein --outfile $2 --outfmt=phy --force
    rm -f RB.seq
}

align_pair () {
    python Exercise04.py -if $1 $2 -of RB.seq
    clustalo -i RB.seq --seqtype=Protein --outfile $3 --outfmt=phy --force
    rm -f RB.seq
}

protdist () {

    phylip protdist <<EOF
$1
2
Y
EOF

    mv outfile $2

}

neighbor () {

    phylip neighbor <<EOF
$1
2
Y
EOF

    mv outfile $2.outfile
    mv outtree $2.outtree

}

Q2 () {
    align TME4bis/RB_sequences/PF01599.12.seq aln-phylip
    protdist aln-phylip protDist_outfile
    neighbor protDist_outfile NJ
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
    neighbor protDist_outfile NJ
    archaeopteryx NJ_outtree
    rm -f aln-phylip
}

Q4 () {
    pfam_path='TME4bis/RB_clades'
    pfam_pair_path='TME4bis/RB_clades_pairs'
    pfam_pair_concat_path='TME4bis/RB_clades_concat_tree'
    pfam_pair_prodist_path='TME4bis/RB_clades_pairs_prodist'
    pfam_pair_tree_path='TME4bis/RB_clades_pairs_tree'

    mkdir -p $pfam_pair_path $pfam_pair_prodist_path $pfam_pair_tree_path $pfam_pair_concat_path

    for pfam in `ls $pfam_path` ; do
	curr_path=`pwd`
	pppath=$pfam_pair_path/$pfam
	mkdir -p $pppath
	cd $pfam_path/$pfam
	clade_list=(`ls *seq`)
	cd $curr_path
	for c1 in `seq 0 $(expr ${#clade_list[@]} - 1)` ; do
	    for c2 in `seq $(expr ${c1} + 1) $(expr ${#clade_list[@]} - 1)` ; do
		clade_file1=$pfam_path/$pfam/${clade_list[$c1]}
		clade_file2=$pfam_path/$pfam/${clade_list[$c2]}
		clade1=(`echo ${clade_list[$c1]} | tr '.' ' '`)
		clade2=(`echo ${clade_list[$c2]} | tr '.' ' '`)
	        align_pair $clade_file1 $clade_file2 $pppath/${clade1[0]}_${clade2[0]}.phylip
	    done
	done
	python concat_phylip.py -if $pppath/*phylip -of $pfam_pair_concat_path/$pfam.phylip && \
	protdist $pfam_pair_concat_path/$pfam.phylip $pfam_pair_prodist_path/$pfam.prodist  && \ 
	neighbor $pfam_pair_prodist_path/$pfam.prodist $pfam_pair_tree_path/$pfam
    done

    python build_input_supertree.py -if $pfam_pair_tree_path/*outtree -of input.trees

    supertree -n output.tree input.trees

    archaeopteryx output.tree
}

Q4_clean () {
    pfam_path='TME4bis/RB_clades'
    pfam_pair_path='TME4bis/RB_clades_pairs/'
    pfam_pair_tree_path='TME4bis/RB_clades_pairs_tree/'
    pfam_pair_concat_path='TME4bis/RB_clades_concat_tree'
    pfam_pair_prodist_path='TME4bis/RB_clades_pairs_prodist/'
    rm -rf $pfam_path $pfam_pair_path  $pfam_pair_prodist_path $pfam_pair_tree_path $pfam_pair_concat_path
    cd TME4bis
    tar -xf RB_clades.tar.gz &> /dev/null 
}

# Q2 
# Q3 
# Q4_clean
# Q4

