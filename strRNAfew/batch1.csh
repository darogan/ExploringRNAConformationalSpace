foreach rna ( RF00028 RF00050 RF00162 RF00167 )
	echo $rna
	cd $rna
	tcsh code/runs.csh 10
	cd ..
end
