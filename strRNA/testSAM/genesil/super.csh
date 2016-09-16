code/sim > sim.log
cat dump.pdb | sort -n -k6 | sed 's/ CA / P  /' > model.pdb
code/scoreWrms pairs.dat model.pdb 2gis.bak | tail -1
rasmol -script chain.ras super.pdb
