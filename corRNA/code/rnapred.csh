# 1 = no of seqs

# make and enter rnapred work sub-directory
if ( -e rnapred ) then
        echo rnapred exists
	cd rnapred
else
        echo rnapred created
        mkdir rnapred
	if ( -e final.seq ) cp final.seq rnapred
	cd rnapred
	ln -s ~/rnapred main
	ln -s ~/util
	ln -s ~/rnacor/corRNA/code
endif
# run RNAfold on each sequence
grep    '>' ../final.afa | head -$argv[1] > rnas.ids
grep -v '>' ../final.afa | head -$argv[1] > rnas.seq
rm pairs.dat
@ n = 0
while ( $n < $argv[1] )
	@ n++
	head -$n rnas.ids | tail -1 
	head -$n rnas.seq | tail -1 > rna.seq 
	echo "@" >> rna.seq
	RNAfold < rna.seq > rna.out
	cat rna.out
	grep '\[[1-9]' rna.ps | grep -v '\.' | tr -d "][" | awk '{print $1,$2}' >> pairs.dat
end
sort pairs.dat | awk -v n=1 '{if($1==a && $2==b){n++;}else{print $1,$2," 0 8 ",0.1*n; a=$1; b=$2; n=1}}' | sort -nr -k 5 > bpairs.dat
tcsh code/view.csh bpairs.dat ../true.pdb 50 5
