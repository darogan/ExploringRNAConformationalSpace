# 1 = RF00123

if ( -e $argv[1] ) then
	if ( -e dump/$argv[1] ) rm -r dump/$argv[1]
	mv hold/$argv[1] dump
	mv $argv[1] hold
endif
mkdir $argv[1]
cp $argv[1].seq $argv[1]/query.seq
cd $argv[1]
ln -s ~/corRNA/code
ln -s ~/corRNA/rfam
ln -s ~/util
ln -s ~/util/data
ln -s ~/util/multal/mulsel
ln -s ~/util/multal/multas
cp rfam/$argv[1].afa.txt .
cp rfam/$argv[1].fa.gz .
gunzip $argv[1].fa.gz
mv $argv[1].fa rfam.fa
grep '>' rfam.fa | wc -l
cat rfam.fa | sed 's/>/+ /' | awk '{if($1=="+"){print "*\n>" $2 "\n" $0}else{print}}' | sed 's/+ //' | sed '1 d' > rfam.seq
echo "*" >> rfam.seq
@ n = `grep '>'  rfam.seq | wc -l`
echo $n rfam sequences
echo Setup done
if ( -e query.seq ) then
	@ u = `tail -2 query.seq | grep U | wc -c`
	if ( $u > 0) then
		echo query.seq contains Us
		exit
	endif
	echo running now
else
	echo query.seq not found
	exit
endif
sleep 3
tcsh code/rerun.rfam.csh 4000
echo Aligning to seed profile
tcsh code/addseed.csh $argv[1] 4
