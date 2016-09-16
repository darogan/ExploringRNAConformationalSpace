# 1 = RF00123

if ( -e $argv[1] ) then
	rm -r $argv[1].hold
	mv $argv[1] $argv[1].hold
endif
rm -r $argv[1]
mkdir $argv[1]
cp $argv[1].seq $argv[1]/query.seq
cd $argv[1]
ln -s .. home
ln -s ~/util
ln -s ~/util/data
ln -s ~/util/multal/mulsel
ln -s ~/util/multal/multas
cp ~/Downloads/$argv[1].fa.gz .
gunzip $argv[1].fa.gz
mv $argv[1].fa rfam.fa
grep '>' rfam.fa | wc -l
cat rfam.fa | sed 's/>/+ /' | awk '{if($1=="+"){print "*\n>" $2 "\n" $0}else{print}}' | sed 's/+ //' | sed '1 d' > rfam.seq
echo "*" >> rfam.seq
@ n = `grep '>'  rfam.seq | wc -l`
echo $n rfam sequences
if ( -e query.seq ) then
	echo running
else
	echo query.seq not found
	exit
endif
sleep 3
tcsh home/rerun.rfam.csh
