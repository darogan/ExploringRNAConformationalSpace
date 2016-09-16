cat RF00017.fa | sed 's;[/-]; ;g' | awk '{if(match($1,">")){n =$3-$2; if(n<0){n=-n}; m=0; print $1,n}else{m++; if(n>200 && m<3){out=0}else{out=1}; if(out){print $1}}}' >chop.fa
@ n = `grep '>' chop.fa | wc -l`
echo $n rfam sequences
if ( $n > 5000 ) then
	echo reducing with mulsel
	cp pdb.seq final.seq
	cat chop.fa | awk '{n++; if(match($1,">")){if(n>1){print "*"}; print $1 "\n" "rfam"}else{print $1}}' >> final.seq
exit
	tcsh code/mulsel.csh 6000 $argv[1]
	echo ">PDBseq" > rfam.fa
	tail -1 pdb.fa >> rfam.fa
	cat final.seq | sed 's/SEED://' | grep -v rfam | grep -v ']' | tr -d "*" >> rfam.fa
endif
