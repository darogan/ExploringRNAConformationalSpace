# 1 = No of models, 2 = max No of pairs

cp code/noview.run test.run

head -$argv[2] pairs.grem.dat > pairs.dat
@ n = `cat pairs.dat | wc -l`
echo Running with $n grem pairs
tcsh code/models.csh $argv[1] > pairs.grem.rms
sleep 1
mv models.pdb models.grem.pdb

echo Unsmoothed
cat pairs.grem.rms | awk '{n++; s+=$3; printf("%7.3f\n",s/n)}' | tail -1

echo Smoothed1x
cat pairs.grem.rms | awk '{n++; s+=$6; printf("%7.3f\n",s/n)}' | tail -1

echo Smoothed2x
cat pairs.grem.rms | awk '{n++; s+=$9; printf("%7.3f\n",s/n)}' | tail -1

code/score models.grem.pdb | grep avd > pairs.grem.plot
