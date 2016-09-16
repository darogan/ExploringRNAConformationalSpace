grep '^[1-9A-Z]' blast162.out | grep -v '                         [1-9]' > blast.out1
cat blast.out1 | sed 's/^/>/' | sed 's/>........  ........ /&#/' | tr "#" "\n" | cut -c 1-94 | sed 's/ /-/g' > blast.out.afa
grep '^>' blast.out.afa | wc -l
more blast.out.afa
