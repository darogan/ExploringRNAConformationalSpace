# 1 = nn (03 = PUZ0003)

# make PUZ00nn.txt and .seq from: http://ahsoka.u-strasbg.fr/rnapuzzles/problems_past.php

rm -r PUZ00$argv[1]

echo Setup work directory
tcsh code/setup.work.csh PUZ00$argv[1]
cd PUZ00$argv[1]
mv ~/Downloads/rfamseq12.chal_$argv[1].aln .
ls -l rfamseq12.chal_$argv[1].aln
sleep 1
echo
echo Convert alignment format
tcsh code/stock2fast.csh rfamseq12.chal_$argv[1].aln
wc rfam.afa
sleep 1
echo
echo Remove sequences with large gaps over 1/4
tcsh code/cutgap.csh rfam 4
wc rfam.cut.afa
sleep 1
echo
echo Run mulsel
tcsh code/makeafa.csh PUZ00$argv[1] 4000
wc final.afa
