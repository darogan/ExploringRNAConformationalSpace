@ n = 0
while ( $n < $argv[1] )
	@ n++
	code/sim full.run > test.log
	grep ' E ' dump.pdb | sort -n -k6 >> dumpE.pdb
	echo TER >> dumpE.pdb
end
rasmol -script chain.ras dumpE.pdb
