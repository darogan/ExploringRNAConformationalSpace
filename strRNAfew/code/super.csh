# 1 = smooth level to keep, 2 = smooth level to view

cat dump.pdb | sort -n -k6 | sed 's/ CA / P  /' | sed 's/XXX/  G/' > model.pdb
code/superWrms pairs.dat model.pdb true.pdb secs.dat | tail -1
cp super$argv[1].pdb super.pdb
grep ' A ' super.pdb > model.pdb
if ( $#argv < 2 ) exit
rasmol -script code/super.ras super$argv[2].pdb &
