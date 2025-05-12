set term @PNG enhanced size 800,400

! awk '{ if ($1 == "file:") file = $2; else print $0 > file; }' < out

unset key
unset xtics
unset ytics
set xrange [-0.5:0.5]
set yrange [0:1.2]
set multiplot layout 3,3 scale 1.05,1.05
plot 'eta-0' u 1:3:($2+$3) w filledcu lc 3, '' u 1:3 w l lc 1
plot 'eta-1' u 1:3:($2+$3) w filledcu lc 3, '' u 1:3 w l lc 1
plot 'eta-2' u 1:3:($2+$3) w filledcu lc 3, '' u 1:3 w l lc 1
plot 'eta-3' u 1:3:($2+$3) w filledcu lc 3, '' u 1:3 w l lc 1
plot 'eta-4' u 1:3:($2+$3) w filledcu lc 3, '' u 1:3 w l lc 1
plot 'eta-5' u 1:3:($2+$3) w filledcu lc 3, '' u 1:3 w l lc 1
plot 'eta-6' u 1:3:($2+$3) w filledcu lc 3, '' u 1:3 w l lc 1
plot 'eta-7' u 1:3:($2+$3) w filledcu lc 3, '' u 1:3 w l lc 1
plot 'eta-8' u 1:3:($2+$3) w filledcu lc 3, '' u 1:3 w l lc 1
unset multiplot

! rm -f eta-? level-?
