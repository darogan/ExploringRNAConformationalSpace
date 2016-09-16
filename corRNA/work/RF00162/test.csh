foreach seq ( `cat mulsel.keep` )
echo $seq
	@ at = $seq * 2
	head -$at rfam.afa | tail -2
end
