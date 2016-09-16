mv ~/Downloads/RF0$argv[1].afa.txt RF0$argv[1].afa.txt
cat RF0$argv[1].afa.txt | grep '^>' | wc -l
cat RF0$argv[1].afa.txt | sed 's/U/T/g' | sed 's/>/#>/' | sed 's/[0-9]$/&#/' | tr -d "\n" | tr "#" "\n" > RF0$argv[1].afa
sed -i '1 d' RF0$argv[1].afa
tcsh cutgap.csh RF0$argv[1] 4
