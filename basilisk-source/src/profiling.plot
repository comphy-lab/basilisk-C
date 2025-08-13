unset mouse
do for [i=0:1000000] {
  load ("< awk -f $BASILISK/trace.awk < profiling")
  pause 10
}
