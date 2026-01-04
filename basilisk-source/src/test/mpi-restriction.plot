set size ratio -1
unset key
unset xtics
unset ytics
plot './out' w l, './out-1' w l, './out-2' w l, './log' u 1:2:3 w labels, './log-1' u 1:2:3 w labels, './log-2' u 1:2:3 w labels 
