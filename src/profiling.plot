unset mouse
load ("< awk -f $BASILISK/trace.awk < profiling")
pause 10
reread
