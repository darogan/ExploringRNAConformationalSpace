# 1 = RF00123 (dir exists and in it)

rm temp*
echo Extracting $argv[1]
cat ../Rfam.full | awk -v s=$argv[1] '{if($2=="AC"){if($3==s){out=1}else{out=0}}; if(out==1){print}}' > $argv[1].sto
grep -v '#' $argv[1].sto | awk '{if($1==""){n++; m=0}; m++; print m "." n, $2}' | sort -n | sed 's/\.1 /&#/' | awk '{print $2}' | tr -d "\n" | tr "#" "\n" > temp1.seq
if ( -e query.gap ) cp query.gap temp2.seq
grep -v '^$' temp1.seq | sed 's/\./-/g' | sed 's/U/T/g' >> temp2.seq
cat temp2.seq | awk '{n++; print ">rfam" n "\n" $1}' > $argv[1].afa
cat temp2.seq | awk -v s=$argv[1] '{n++; print ">RFAM" n "\n" s,n "\n" $1 "*"}' > $argv[1].ali
cat $argv[1].ali | tr -d "-" > $argv[1].seq
cat temp2.seq | tr "-" "x" | awk '{n++; print n "\t" $1}' > $argv[1].xgaps
cat temp2.seq | tr -d "-"  | awk '{n++; print n "\t" $1}' > $argv[1].nogap
echo $argv[1].xgaps
head $argv[1].xgaps
echo $argv[1].nogap
head $argv[1].nogap
echo $argv[1].afa
head $argv[1].afa
echo $argv[1].ali
head $argv[1].ali
echo $argv[1].seq
head $argv[1].seq
rm temp*
