set size ratio -1
unset key
plot 'out' w l, 'log' u 1:2:(5-$3) pt 6 lc 3 ps variable
