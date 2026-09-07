#!/bin/bash
# Finish the threshold-23 tower: level 12->13 output is already being produced in
# t23run/.  This waits for it, dedups, then climbs 13->14 and 14->15.
set -e
cd "$(dirname "$0")"
export OMP_NUM_THREADS=1

while [ "$(grep -l raw_extensions_kept t23run/log*.txt 2>/dev/null | wc -l | tr -d ' ')" -lt 16 ]; do
  sleep 15
done
echo "LEVEL 12->13 shards complete"
grep -h raw_extensions_kept t23run/log*.txt | awk '{for(i=1;i<=NF;i++){if($i~/^graphs=/){split($i,a,"=");g+=a[2]};if($i~/^raw_extensions_kept=/){split($i,b,"=");h+=b[2]}}} END{print "  sources",g,"raw_extensions",h}'
cat t23run/r*.g6 > T13_raw22.g6
labelg -q T13_raw22.g6 | sort -u > T13_thresh22.g6
echo "  |T13(M<=22)| = $(wc -l < T13_thresh22.g6)  (raw $(wc -l < T13_raw22.g6))"

mkdir -p t23run14
find t23run14 -type f -delete
seq 0 15 | xargs -P 16 -I{} sh -c 'python3 expand_prune.py T13_thresh22.g6 12 23 t23run14/r$1.g6 --res $1 --mod 16 --K 128 > t23run14/log$1.txt 2>&1' _ {}
echo "LEVEL 13->14 complete"
cat t23run14/r*.g6 > T14_raw23.g6
labelg -q T14_raw23.g6 | sort -u > T14_thresh23.g6
echo "  |T14(M<=23)| = $(wc -l < T14_thresh23.g6)  (raw $(wc -l < T14_raw23.g6))"

mkdir -p t23run15
find t23run15 -type f -delete
seq 0 15 | xargs -P 16 -I{} sh -c 'python3 expand_prune.py T14_thresh23.g6 13 23 t23run15/r$1.g6 --res $1 --mod 16 --K 128 > t23run15/log$1.txt 2>&1' _ {}
echo "LEVEL 14->15 complete"
cat t23run15/r*.g6 > T15_raw23.g6
echo "  order-15 matrices with M <= 23, raw count = $(wc -l < T15_raw23.g6)"
if [ -s T15_raw23.g6 ]; then
  labelg -q T15_raw23.g6 | sort -u > T15_thresh23.g6
  echo "  |T15(M<=23)| = $(wc -l < T15_thresh23.g6)  => F(15) = 23"
else
  echo "  T15(M<=23) EMPTY => F(15) >= 25"
fi
echo TOWER23_DONE
