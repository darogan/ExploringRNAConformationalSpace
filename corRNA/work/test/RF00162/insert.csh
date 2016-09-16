# insert query seqs  (every 20)
#cat mulsel.hold | sed 's/>/> /' | awk '{if($1==">"){n++; if(n==20){print "here"; n=0} print ">" $2; }else{ print }}' > temp.seq
cat seed.seq | sed 's/>/> /' | awk '{if($1==">"){n++; if(n==20){print "here"; n=0} print ">" $2; }else{ print }}' > temp.seq
cat temp.seq | sed 's/^/echo "/' | sed 's/$/"/'  | sed 's/echo "here"/cat insert.seq/' > temp.csh
sed 's/SEED://' query.seq > insert.seq
cp query.seq mulsel.seq
tcsh temp.csh >> mulsel.seq
cat insert.seq >> mulsel.seq
# align
./multas > multas.log
# strip out queries
head -2 final.afa > strip.afa
awk '{if(match($1,">")){  if(match($0,"RF00162")){skip=1}else{skip=0} }; if(skip==0){print}}' final.afa >> strip.afa
mv strip.afa final.afa
# count
grep seqs final.aln
grep '>' final.afa | wc -l
echo top
grep -v '>' final.afa | head
echo bot
grep -v '>' final.afa | tail
