#!/bin/tcsh
#
# 1 = list of pairs, 2 = PDB file, 3 = number of pairs to view, [4 = min seq. separation]
#
if ( $#argv > 3 ) then
	echo Excluding separations under $argv[4]
endif
# basepairs
if ( $#argv > 3 ) then
	awk -v n=$argv[4] '{if($2-$1 > n+1){ print}}' rnapred/bpairs.dat | head -$argv[3] > link.tmp
else
	head -$argv[3] rnapred/bpairs.dat > link.tmp
endif
cat link.tmp | awk '{printf("CONECT %4d    0    0    0    0 %4d\n", $1, $2)}' > link.dat
# contacts
if ( $#argv > 3 ) then
	awk -v n=$argv[4] '{if($2-$1 > n+1){ print}}' $argv[1] | head -$argv[3] > bond.tmp
else
	head -$argv[3] $argv[1] > bond.tmp
endif
# mark common pairs
cat link.tmp | awk '{print $0,$1*$2}' >  both.tmp
cat bond.tmp | awk '{print $0,$1*$2}' >> both.tmp
sort -n -k6 both.tmp > sort.tmp
cat sort.tmp | awk '{if($6==last){print $0,"M"}else{print $0, "N"}; last=$6}' | grep '0 8' > bond.tmp
# resort on score and split
cat bond.tmp | sort -nr -k5 | awk '{n++; print $0, n}' > bond0.tmp
grep 'M' bond0.tmp > bondM.tmp
grep 'N' bond0.tmp > bondN.tmp
cat bondM.tmp | awk '{m=100-$8*4; if(m>0){print $1, $2, m}}' > bond1.tmp
cat bond1.tmp | awk '{printf("select %d,%d\nhbond %d\ncolour hbonds blue \n", $1, $2, $3)}' > bond1.ras
cat bondN.tmp | awk '{m=100-$8*4; if(m>0){print $1, $2, m}}' > bond2.tmp
cat bond2.tmp | awk '{printf("select %d,%d\nhbond %d\ncolour hbonds green\n", $1, $2, $3)}' > bond2.ras
#
cat bond[12].tmp | awk '{printf("CONECT %4d    0    0    0    0 %4d\n", $1, $2)}' >> link.dat
cat $argv[2] link.dat > link.pdb
cat viewRNA.ras bond1.ras bond2.ras > temp.ras
echo select >> temp.ras
rasmol -script temp.ras link.pdb
rm *.tmp
