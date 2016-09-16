# 1 = file, 2 = chain ID (otherwise all)

cat $argv[1].pdb | grep -v REMARK | grep ' P  '  > temp
awk '/ATOM/ {n++; printf("ATOM %6d  P     %s A%4d     %7.3f %7.3f %7.3f %5.2f %5.2f\n", n,$4,n,$7,$8,$9,$10,$11)}' temp > temp.cas
if ( $#argv > 1 ) then
	eval "grep ' $argv[2] ' temp.cas" > $argv[1].cas
	rm temp.cas
else
	echo no chain ID
	mv temp.cas $argv[1].cas
endif
