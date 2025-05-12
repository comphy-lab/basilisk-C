! awk '{ if ($2 == "file:") file = $3; else print $0 > file; }' < out
set term @PNG enhanced size 640,426

set size ratio -1
unset key
unset xtics
unset ytics
unset border
unset colorbox
set pm3d
set pm3d map interpolate 1,1
set xrange [-500:500]
set yrange [-500:500]
set zrange [-0.05:0.05]
set palette rgbformulae 33,13,10

set multiplot layout 2,3 scale 1.6,1.6
splot 't-0'
splot 't-300'
splot 't-600'
splot 't-900'
splot 't-1200'
splot 't-1500'
unset multiplot

reset
set term @PNG enhanced size 426,520
set output 'diff.png'

set pm3d map
set palette rgbformulae 33,13,10
set xrange [-500:500]
set yrange [-500:500]
set zrange [-0.0002:0.0002]
splot 'diff' using 1:2:7 with pm3d
 
! rm -f t-* level-*
