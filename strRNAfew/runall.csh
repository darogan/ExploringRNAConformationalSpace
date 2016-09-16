foreach rna ( `ls -d RF*` )
	echo $rna
	cd $rna
	tcsh code/runs.csh 10
	cd ..
end
