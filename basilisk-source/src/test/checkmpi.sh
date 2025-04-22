# checks the consistency of send/receive buffers

echo "checkmpi" > /dev/stderr

t=""
npe=`ls mpi-level-rcv-$t* | wc -l | awk '{print $1 - 1}'`
for i in `seq 0 1 $npe`; do
    for j in `seq 0 1 $npe`; do
	for op in level level-root restriction; do
	    awk -v j=$j '{if ($4 == j) print $1,$2,$3;}' \
		< mpi-$op-rcv-$t$i > rcv-$op-$i-$j
	    awk -v i=$i '{if ($4 == i) print $1,$2,$3;}' \
		< mpi-$op-snd-$t$j > snd-$op-$j-$i
	    if ! diff rcv-$op-$i-$j snd-$op-$j-$i > diff; then
		echo \'rcv-$op-$i-$j\', \'snd-$op-$j-$i\' > /dev/stderr
		cat diff > /dev/stderr
	    else
		rm -f rcv-$op-$i-$j snd-$op-$j-$i
	    fi
	    rm -f diff
	done
    done
done

cat depth-$t* | sort -r | awk '
/^depth:/ { depth[$2] = $3; }
{
  if ($1 != "=======")
    if ($3 > depth[$2])
      printf ("process %d thinks process %d is deeper than it really is! (%d > %d)\n",
              $1, $2, $3, depth[$2]) > "/dev/stderr";
}
'
