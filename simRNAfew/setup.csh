# 1 = RF00123/1xyz.pdb

rm code
ln -s ../code
ln -s code/data
ln -s .. home
cp ../work/$argv[1] .
cp ../work/$argv[1] full.pdb
head full.pdb
echo > bracket.dat
echo run GeneSilico analyse with renum and secstr
# http://iimcb.genesilico.pl/modernaserver/submit/analyse/
echo mv ~/Downloads/123.cleaned.renumbered_chain.pdb full.pdb
