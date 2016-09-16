foreach rna ( `ls -d RF* | grep -v ct` )
	echo
	echo $rna
	code/scoreWrms $rna/pairs.true.dat $rna/true.pdb
	code/scoreWrms $rna/pairs.pred.dat $rna/true.pdb
	@ len = `cat $rna/true.pdb | wc -l`
	set grem = `code/scoreWrms $rna/gremlin.300 $rna/true.pdb | tail -1`
	echo $grem   len = $len
end
