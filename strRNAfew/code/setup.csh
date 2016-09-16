# 1 = RF00123

# run genesilico
#
# http://iimcb.genesilico.pl/modernaserver/submit/analyse/
#
# load full pdb and click "get secondary structure"
# wait
# click on the applet view box
# select save as CT format
# (and eps)
#

if ( -e $argv[1].ct ) then
	echo Setting-up dir $argv[1]
else
	echo Need to run Genesilico to get basepairs
	exit
endif

mkdir $argv[1]
cp $argv[1].ct $argv[1]/pairs.ct
cd $argv[1]
ln -s ~/rnacor/strRNA/code

# write basebairs
code/basepairs
#gnuplot> plot [1:94][1:94] 'pairs.out' not, x not, 'secs.out' u 1:1 not

# find stems and write in segs.dat
@ len = `head -1 pairs.ct`
@ gap = 5 + $len / 100
echo gap penalty = $gap
tcsh code/boxplot.csh 1 $gap $len
cp line.plot segs.dat
cp line.plot segs.raw
# probably need to edit ends

# write domains for sim in stems.dat (stems.pdb)
cp ~/rnacor/corRNA/work/$argv[1]/true.pdb .
code/stems pairs.dat true.pdb
rasmol -script code/chain.ras stems.pdb &
cp stems.dat stems.hold
cat segs.dat

cat pairs.dat | sed 's/ 0 8  1.0/   0  10000/' > pairs.true.dat
cat  ~/rnacor/corRNA/work/$argv[1]/rnapred/bpairs.dat | sed 's/ 0 8  0.1/   0  10000/' > pairs.pred.dat
head -100 ~/rnacor/corRNA/work/$argv[1]/gremlin.300 | sed 's/ 0 8  0./   0  /' | sed 's/.$//' > pairs.grem.dat
# pairs.dat read in test.run
cp pairs.true.dat pairs.dat

cp code/nomove.run test.run

exit

# check length and look for diff
grep ATOM stems.dat | wc -l
grep ATOM stems.dat | sort -n -k6 | awk '{ if($6==last){ print "********"} print $0; last=$6}' | more

echo TER >> stems.dat
awk '{print "REBOND pdbid",$1,$1+1}' breaks.dat >> stems.dat
echo RETERM pdbid 1 $len >> stems.dat

# use +NORUN +NOMOVE first then -NORUN +NOMOVE for short test run
code/sim
# then check superposition
tcsh code/super.csh 0 0

# run with pairs.true.dat and move = 0 --> new true.pdb

tcsh code/pairs.csh 20

tcsh runs.csh

gnuplot> plot 'save10/pairs.true.plot' u 3:9 w l, 'save10/pairs.pred.plot' u 3:9 w l, 'save10/pairs.grem.plot' u 3:9 w l
gnuplot> plot 'save10/pairs.true.plot' u 3:6 w l, 'save10/pairs.pred.plot' u 3:6 w l, 'save10/pairs.grem.plot' u 3:6 w l

