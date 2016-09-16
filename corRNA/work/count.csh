foreach rna ( `cat used.list` )
	@ leng = `cat $rna/true.pdb | wc -l`
	@ full = `grep '>' ../infernal/$rna/rfam.fa | wc -l`
	@ muls = `grep '>' $rna/rfam.afa | wc -l`
	@ grem = `cat $rna/gremlin.300 |  wc -l`
	@ secs = `grep bonds $rna/mark.ras | wc -l`
	@ gren = `grep green $rna/mark.ras | wc -l`
	@ blue = `grep blue  $rna/mark.ras | wc -l`
	@ reds = `cat $rna/mark.ras | grep hbonds | grep red | wc -l`
	@ yelo = `cat $rna/mark.ras | grep hbonds | grep yellow | wc -l`
	@ link = `cat $rna/rnapred/bpairs.dat | wc -l`
	@ true = `cat  ~/rnacor/strRNA/$rna/pairs.true.dat | wc -l`
	@ pred = `cat  ~/rnacor/strRNA/$rna/pairs.pred.dat | wc -l`
	@ cut9 = `cat $rna/gremlin.raw | grep -v prob | awk '{if($7>0.9){n++}; print n}' | tail -1`
	@ cut7 = `cat $rna/gremlin.raw | grep -v prob | awk '{if($7>0.7){n++}; print n}' | tail -1`
	@ cut5 = `cat $rna/gremlin.raw | grep -v prob | awk '{if($7>0.5){n++}; print n}' | tail -1`
	echo $rna $leng $full $muls $grem $secs $gren $blue $reds $yelo $link $true $pred $cut9 $cut7 $cut5
end
