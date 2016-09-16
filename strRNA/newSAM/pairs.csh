cp noview.run test.run
echo basepairs
cp bpairs.dat pairs.dat
tcsh models.csh $argv[1] > bpairs.rms
mv models.pdb models.sec.pdb
echo corpairs
cp cpairs.dat pairs.dat
tcsh models.csh $argv[1] > cpairs.rms
mv models.pdb models.cor.pdb
