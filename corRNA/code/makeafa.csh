# 1 = RF00123, 2 = min kept (4000)

if ( -e home/infernal/$argv[1]/rfam.cut.afa ) then
	cp home/infernal/$argv[1]/rfam.cut.afa rfam.tmp
else
	cp rfam.cut.afa rfam.tmp
endif
cat rfam.tmp | grep -v seq | tr "U" "T" | awk '{n++; print ">seq" n "\n" $1}' > rfam.afa
# set up mulsel run
cat rfam.afa | awk '{print $1; if(match($1,">")){print "rfam"}else{ print "*"}}' | sed 's/>seq1$/SEED:>seq1/' > final.seq
tcsh code/mulsel.csh $argv[2] 4
# compile reduced afa
grep seq final.seq | sed 's/>seq//' > mulsel.keep
head -2 rfam.afa > mulsel.afa
foreach seq ( `cat mulsel.keep | grep -v SEED` )
	@ at = $seq * 2
	head -$at rfam.afa | tail -2 >> mulsel.afa
end
cp mulsel.afa final.afa
@ n = `grep '>' final.afa | wc -l`
echo $n sequences in final.afa
head final.afa
