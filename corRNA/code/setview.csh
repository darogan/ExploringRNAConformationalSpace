if ( ! -e gremlin.raw ) then
	echo Paste gremlin results
	cat > gremlin.raw
endif
tcsh code/gremlin2psicov.csh
if ( -e rnapred ) then 
	echo rnapred exists
else
	tcsh code/rnapred.csh 1
endif
tcsh code/viewRNA.csh gremlin.300 true.pdb 100 5
