# 1 = gapped sequences, 2 = fractional gap cutoff

rm $argv[1].cut.afa
rm cols.dat
grep '^>' $argv[1].afa > codes.dat
echo `cat codes.dat | wc -l` starting sequences

set len = `cat $argv[1].afa | head -2 | tail -1 | wc -c`
echo len = $len

head -2 $argv[1].afa | tail -1 | tr "-" "x" > query.seq
touch cols.dat
@ n = 0
while ( $n < $len )
	@ n++
	# test query sequence for gap (assume first)
	set res = `cut -c $n query.seq`
	if ( $res == "x" ) continue
	# if not gap, cut out and add column
	grep -v '>' $argv[1].afa | cut -c $n > col.tmp
	paste -d '' cols.dat col.tmp > cols.tmp
	mv cols.tmp cols.dat
end
paste -d '#' codes.dat cols.dat | tr " " "~" > seqs.tmp
set len = `cat seqs.tmp | head -2 | tail -1 | wc -c`
echo Query length is $len
@ m = 1
@ max = $len / $argv[2]
echo Removing sequences with over $max gaps
foreach seq (`cat seqs.tmp`)
	set n = `echo "$seq" | tr -d "-" | wc -c`
	@ gaps = $len - $n
	if ( $gaps > $max && $m > 1 ) continue
	echo "$seq" | tr "#" "\n" | tr "~" " " >> $argv[1].cut.afa
	@ m++
end
echo $m left
