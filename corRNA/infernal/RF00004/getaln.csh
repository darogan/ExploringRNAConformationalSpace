#cmalign rfam.cm rfam.fa > temp.out
cat temp.out | grep -v '=G[SRFC]' | tr "a-z." "A-Z-" | sed '1,2 d' | sed '$ d' > temp.txt
cat temp.txt | awk '{if($1==""){n++; m=0}; m++; print m "." n, $2}' | sort -n | sed 's/\.1 /&#/' | awk '{print $2}' | tr -d "\n" | tr "#" "\n" > temp.raw
cat temp.raw | grep -v '^$' | awk '{n++; print ">seq" n "\n" $1}' > rfam.afa
tcsh ../cutgap.csh rfam 4
rm temp.*
