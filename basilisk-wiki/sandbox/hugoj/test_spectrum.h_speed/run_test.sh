#!/usr/bin/env bash

# Monitor with htop and nvtop

set -x
SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &>/dev/null && pwd)

echo $SCRIPT_DIR
NAME=speed_wave_ini

make -k clean
rm -r $NAME

list_N=(128 256 512)
list_nl=(10 30 50 70 90)
NLD=5
ND=64
# run once for compile and default res
CFLAGS="-DTRACE=2" make $NAME.tst
mv $NAME/out $NAME/out_N${ND}_nl${NLD}
mv $NAME/log $NAME/log_N${ND}_nl${NLD}

cd $NAME

# run over N
for N in "${list_N[@]}"; do
  NL=$NLD
  ./$NAME $N $NLD 2>log >out
  mv out out_N${N}_nl${NL}
  mv log log_N${N}_nl${NL}
done

# run over nl
for NL in "${list_nl[@]}"; do
  N=$ND
  ./$NAME $N $NL 2>log >out
  mv out out_N${N}_nl${NL}
  mv log log_N${N}_nl${NL}
done
