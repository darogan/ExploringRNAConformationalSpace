#!/bin/tcsh
#
# 1 = list of pairs, 2 = PDB file, 3 = number of pairs to view, [4 = min seq. separation]
#
if ( $#argv > 3 ) then
	echo Excluding separations under $argv[4]
	awk -v n=$argv[4] '{if($2-$1 > n+1){ print}}' $argv[1] | head -$argv[3] > temp1.dat
else
	head -$argv[3] $argv[1] > temp1.dat
endif
cat temp1.dat | awk '{printf("CONECT %4d    0    0    0    0 %4d\n", $1, $2)}' > temp2.dat
cat temp1.dat | awk '{n++; m=100-n*4; if(m>0){printf("select %d,%d\nhbond %d\n", $1, $2, m)}}' > bond.ras
cat $argv[2] temp2.dat > link.pdb
rm temp[12].dat
cat home/view.ras bond.ras > temp.ras
echo select >> temp.ras
rasmol -script temp.ras link.pdb &
