unset xtics
unset ytics
unset key
unset border
set size ratio -1
plot 'out' w l, 'log' u 1:2
