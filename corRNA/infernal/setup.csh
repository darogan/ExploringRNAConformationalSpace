# 1 = RF00123, 2 = mulsel range  (+|- 1/x)

if ( $#argv < 3 ) then
	@ keeps = 6000
else
	@ keeps = $argv[3]
endif
echo mulsel target = $keeps

if ( $#argv < 2 ) then
	@ range = 4
else
	@ range = $argv[2]
endif
echo mulsel range = plus/minus length/$range

mkdir $argv[1]
cd $argv[1]
ln -s ~/corRNA home
ln -s home/code
cp home/work/$argv[1].seq pdb.seq
# pdb seq top
echo ">"$argv[1] > pdb.fa
tail -2 pdb.seq | head -1 >> pdb.fa
cp pdb.fa rfam.fa
# add rest
cp home/rfam/$argv[1].fa.gz .
gunzip $argv[1].fa.gz
@ n = `grep '>' $argv[1].fa | wc -l`
echo $n rfam sequences
if ( $n > $keeps ) then
	ln -s ~/util/data
	echo reducing with mulsel
	cp pdb.seq final.seq
	cat $argv[1].fa | awk '{n++; if(match($1,">")){if(n>1){print "*"}; print $1 "\n" "rfam"}else{print $1}}' >> final.seq
	tcsh code/mulsel.csh $keeps $range
	echo ">PDBseq" > rfam.fa
	tail -1 pdb.fa >> rfam.fa
	cat final.seq | sed 's/SEED://' | grep -v rfam | grep -v ']' | tr -d "*" >> rfam.fa
else
	cat $argv[1].fa >> rfam.fa
endif
#
# extract cm
cp ../head.cm rfam.cm
cat ../Rfam.cm | awk -v n=$argv[1] '{if(match($2,n)){out++}; if(match($1,"//") && out==2){out=0}; if(out){print}}' >> rfam.cm
echo "//" >> rfam.cm
