foreach rna ( RF00174 RF00380 RF01051 RF01831 )
	echo $rna
	cd $rna
	tcsh code/runs.csh 10
	cd ..
end
