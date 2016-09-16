# reformat gremlin pairs to look like psicov (excluding i-j<5)
# scores are rescaled by rank (using protein length)

@ len = `cat pdb.seq | sed '1,2 d' | sed 's/\*//' | grep -v '^$' | wc -c`
@ len--
echo $len = protein length
grep -v prob gremlin.raw | awk '{n=5; if($2-$1 > n-1){print  $1,$2," 0 8 ",$7}}' > gremlin.fix
cat gremlin.fix | awk -v n=$len '{x++; e = exp(-x*x*0.01/n); printf("%4d %4d  %d %d  %7.5f\n", $1,$2,$3,$4,e)}' > gremlin.dat
head -300 gremlin.dat > gremlin.300
