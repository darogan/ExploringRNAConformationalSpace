foreach rna ( RF00010 RF00059 RF00234 RF00168 )
	echo $rna
	cd $rna
	tcsh code/runs.csh 10
	cd ..
end
