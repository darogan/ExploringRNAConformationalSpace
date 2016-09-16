cmalign rfam.cm rfam.fa > temp.out
cat temp.out | grep -v '=G[SRFC]' | tr "a-z." "A-Z-" | grep -v STOCKHOLM | grep -v WARNING | sed '$ d' > temp.txt
head -2 temp.txt | tail -1 | grep '^$'
if ( $? == 0 ) then
	echo double blank fixed
	sed -i '1 d' temp.txt
endif
cat temp.txt | awk '{if($1==""){n++; m=0}; m++; print m "." n, $2}' | sort -n | sed 's/\.1 /&#/' | awk '{print $2}' | tr -d "\n" | tr "#" "\n" > temp.raw
cat temp.raw | grep -v '^$' | awk '{n++; print ">seq" n "\n" $1}' > rfam.afa
tcsh ../cutgap.csh rfam 4
