awk '/^ATOM / {n++; printf("ATOM %6d  P     %s A%4d     %7.3f %7.3f %7.3f %5.2f %5.2f\n", n,$4,n,$7,$8,$9,$10,$11)}' $argv[1] > $argv[2]
