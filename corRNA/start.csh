# 1 = RF00123, 2 = pdb code

if ( -e pdb/$argv[2].pdb ) then
	echo pdb $argv[2] found
else
	echo no pdb
	exit
endif

if ( -e work/$argv[1].seq ) then
	echo sequence found
else
	echo "SEED:>"$argv[1] > work/$argv[1].seq
	eval "grep $argv[2] seqs/*" | tr "#" "\n" | sed 's/U/T/g' >>  work/$argv[1].seq
	vim  work/$argv[1].seq
endif

if ( -e rfam/$argv[1].fa.gz ) then
	echo Rfam fasta sequences found
else
	if ( -e ~/Downloads/$argv[1].fa.gz ) then
		mv ~/Downloads/$argv[1].fa.gz rfam
		echo copied from Downloads
	else
		echo no Rfam fasta sequences
		exit
	endif
	rm Downloads
endif

echo all setup

cd infernal
tcsh setup.csh $argv[1]
echo check OK
sleep 3
cd $argv[1]
tcsh ../makealn.csh
echo check OK
sleep 3

cd ~/corRNA/work
tcsh code/setup.work.csh $argv[1]
cd $argv[1]
tcsh code/makeafa.csh $argv[1] 2000
head final.afa

cp home/pdb/$argv[2].pdb .
cat pdb.seq

exit
# submit gremlin job
# fix pdb

tcsh code/pdb2pho.csh $argv[2]
tcsh code/setview.csh 
