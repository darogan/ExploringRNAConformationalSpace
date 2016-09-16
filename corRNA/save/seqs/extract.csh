cat pdb_seqres.txt | sed 's/$/ #/' | tr -d "\n" | tr ">" "\n" | grep 'UUU' | grep 'mol:na' | sed 's/:/ /g' > pdbrna.seq
awk '{if($5<500 && $5>50){print}}' pdbrna.seq > pdbrna.50-500.one
awk '{print ">" $1 " #" $0 "*"}' pdbrna.50-500.one | sed 's/ #/#/g' | tr "#" "\n" > pdbrna.50-500.seq

