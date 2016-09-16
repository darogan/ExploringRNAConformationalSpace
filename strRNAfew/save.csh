foreach rna ( `ls -d RF*` )
	echo $rna
	cd $rna
	mkdir lineNruns
	mv save* lineNruns
	mv runs.log lineNruns
	cd ..
end
