rm models.pdb
@ n = 0
while ( $n < $argv[1] )
	@ n++
	tcsh super.csh
	grep ' A ' super.pdb >> models.pdb
	echo TER >> models.pdb
end
cat super.pdb >> models.pdb
rasmol -script models.ras models.pdb &
