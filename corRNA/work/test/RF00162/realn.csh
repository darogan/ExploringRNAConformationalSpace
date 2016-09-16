./multas > multas.log
grep seqs final.aln
grep '>' final.afa | wc -l
echo top
grep -v '>' final.afa | head
echo bot
grep -v '>' final.afa | tail

