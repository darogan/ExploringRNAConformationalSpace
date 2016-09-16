# 1 = RF00123

if ( -e $argv[1] ) then
	if ( -e dump/$argv[1] ) rm -r dump/$argv[1]
	mv hold/$argv[1] dump
	mv $argv[1] hold
endif
mkdir $argv[1]
cp $argv[1].seq $argv[1]/query.seq
cd $argv[1]
ln -s ~/rnacor/corRNA home
ln -s ~/rnacor/corRNA/code
ln -s ~/rnacor/corRNA/rfam
ln -s ~/util
ln -s ~/util/data
cp query.seq pdb.seq
