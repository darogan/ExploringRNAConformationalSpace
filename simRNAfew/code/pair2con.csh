awk '{printf("SLOPE A/%d/MB A/%d/MB 5.0 10.0 1.0\nWELL A/%d/MB A/%d/MB 5.0 10.0 1.0\n", $1,$2,$1,$2)}' pairs.dat > pairs.con
