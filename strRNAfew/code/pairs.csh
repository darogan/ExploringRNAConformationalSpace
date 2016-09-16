# 1 = No of models, 2 = max No of pairs

cp code/noview.run test.run

# run true pairs (genesilico)
head -$argv[2] pairs.true.dat > pairs.dat
@ n = `cat pairs.dat | wc -l`
echo Running with $n true pairs
tcsh code/models.csh $argv[1] > pairs.true.rms
sleep 1
mv models.pdb models.true.pdb

# run predicted pairs (RNAfold)
head -$argv[2] pairs.pred.dat > pairs.dat
@ n = `cat pairs.dat | wc -l`
echo Running with $n pred pairs
tcsh code/models.csh $argv[1] > pairs.pred.rms
sleep 1
mv models.pdb models.pred.pdb

#run correlated pairs (gremlin)
head -$argv[2] pairs.grem.dat > pairs.dat
@ n = `cat pairs.dat | wc -l`
echo Running with $n grem pairs
tcsh code/models.csh $argv[1] > pairs.grem.rms
sleep 1
mv models.pdb models.grem.pdb

echo Unsmoothed
cat pairs.true.rms | awk '/RMS/{n++; s+=$3; printf("%7.3f\n",s/n)}' | tail -1
cat pairs.pred.rms | awk '/RMS/{n++; s+=$3; printf("%7.3f\n",s/n)}' | tail -1
cat pairs.grem.rms | awk '/RMS/{n++; s+=$3; printf("%7.3f\n",s/n)}' | tail -1

echo Smoothed1x
cat pairs.true.rms | awk '/RMS/{n++; s+=$6; printf("%7.3f\n",s/n)}' | tail -1
cat pairs.pred.rms | awk '/RMS/{n++; s+=$6; printf("%7.3f\n",s/n)}' | tail -1
cat pairs.grem.rms | awk '/RMS/{n++; s+=$6; printf("%7.3f\n",s/n)}' | tail -1

echo Smoothed2x
cat pairs.true.rms | awk '/RMS/{n++; s+=$9; printf("%7.3f\n",s/n)}' | tail -1
cat pairs.pred.rms | awk '/RMS/{n++; s+=$9; printf("%7.3f\n",s/n)}' | tail -1
cat pairs.grem.rms | awk '/RMS/{n++; s+=$9; printf("%7.3f\n",s/n)}' | tail -1

code/score models.true.pdb | grep avd > pairs.true.plot
code/score models.pred.pdb | grep avd > pairs.pred.plot
code/score models.grem.pdb | grep avd > pairs.grem.plot
