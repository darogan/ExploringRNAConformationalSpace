code/sim > sim.log
cat dump.pdb | sort -n -k6 | sed 's/ CA / P  /' > model.pdb
code/scoreWrms pairs.dat model.pdb true.pdb | tail -1
#rasmol -script chain.ras super.pdb &
