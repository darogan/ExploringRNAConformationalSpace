 cp rnadb.pdbs.txt rnadb.pdbs.tmp
 vim rnadb.pdbs.tmp
 cat rnadb.pdbs.tmp | tr -d "\n" | tr "#" "\n" > rnadb.pdbs.one
 vim rnadb.pdbs.one
 cat rnadb.pdbs.one | tr "Å" "\n" > rnadb.pdbs.two
 awk '/[0-9]/ {if($1=="#"){print $0}else{print $1; if($1<50 || $1>500){print " X"}}}' rnadb.pdbs.two | tr -d "\n" | tr "#" "\n" > rnadb.pdbs.list
 grep -v 'X$' rnadb.pdbs.list > rnadb.50-500.list
