# 1 = smoothing cycles (10)

code/smooth true.pdb gremlin.300 $argv[1]
grep CON link.pdb >> smooth.pdb
rasmol -script mark.ras link.pdb &
rasmol -script mark.ras smooth.pdb &
