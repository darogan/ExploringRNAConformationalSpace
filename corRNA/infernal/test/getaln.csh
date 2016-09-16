cat test.out | grep -v '=G[SRFC]' | tr "a-z." "A-Z-" | sed '1,2 d' | sed '$ d' > test.txt
cat test.txt | awk '{if($1==""){n++; m=0}; m++; print m "." n, $2}' | sort -n | sed 's/\.1 /&#/' | awk '{print $2}' | tr -d "\n" | tr "#" "\n" > test.raw
cat test.raw | grep -v '^$' | awk '{n++; print ">seq" n "\n" $1}' > test.afa
