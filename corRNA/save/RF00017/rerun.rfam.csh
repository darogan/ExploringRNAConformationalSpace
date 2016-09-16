cp query.seq final.seq
cat rfam.seq | sed 's/U/T/g' >> final.seq
@ n = `grep '>' final.seq | wc -l`
@ leng = `grep -v '[0-9\*]' query.seq | wc -c`
@ leng--
@ le = $leng / 2
@ ng = $leng + 100 
echo $n starting sequences
echo length = $leng
echo range = $le to $ng
cp home/skip99.run skip.run
eval "sed -i 's/50 500/$le $ng/' skip.run"
head -1 skip.run
echo running mulsel
#./mulsel > mulsel.log
#tail -1 mulsel.log
cp start99.seq final.seq
#reduce
foreach to ( 95 90 80 70 60 50 )
	echo reduce to $to
	cp home/skip$to.run skip.run
	./mulsel > mulsel.log
	tail -1 mulsel.log
	cp final.seq final$to.seq
	@ n = `grep '>' final.seq | wc -l`
	if ( $n < 6000 ) break
end
if ( $n > 3000 ) then
	./mulsel > mulsel.log
	tail -1 mulsel.log
	@ n = `grep '>' final.seq | wc -l`
endif
if ( $n > 3000 ) then
	./mulsel > mulsel.log
	tail -1 mulsel.log
endif
# align
echo aligning
cp query.seq mulsel.seq
grep -v 'SEED:' final.seq >> mulsel.seq
cp home/test.run .
./multas > multas.log
grep seqs final.aln
grep '>' final.afa | wc -l
head final.afa
tail final.afa
