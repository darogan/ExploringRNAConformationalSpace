if ( -e models.pdb ) rm models.pdb
@ n = 0
while ( $n < $argv[1] )
	@ n++
	code/sim test.run 1 > sim.log
	tcsh code/super.csh 1
	grep ' A ' super.pdb >> models.pdb
	echo TER >> models.pdb
end
# add true.pdb and view 
grep ' B ' super.pdb >> models.pdb
rasmol -script code/models.ras models.pdb &
