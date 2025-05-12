! awk '{ if ($1 == "file:") file = $2; else print $0 > file; }' < out
 
set term @PNG enhanced size 640,426

set size ratio -1
unset key
unset xtics
unset ytics
unset border
unset colorbox
set pm3d
set pm3d map interpolate 1,1
set palette defined ( 0 0 0 0.5647, 0.125 0 0.05882 1, 0.25 0 0.5647 1, \
    	              0.375 0.05882 1 0.9333, 0.5 0.5647 1 0.4392, \
		      0.625 1 0.9333 0, 0.75 1 0.4392 0, \
		      0.875 0.9333 0 0, 1 0.498 0 0 )

set multiplot layout 2,3 scale 1.6,1.6
set cbrange [0:1]
splot 'field-0' u 1:2:4
splot 'field-1' u 1:2:4
splot 'field-2' u 1:2:4
splot 'field-3' u 1:2:4
splot 'field-4' u 1:2:4
splot 'field-5' u 1:2:4
unset multiplot

reset
set term @PNG enhanced size 640,480
set output 'rate.png'
set xlabel 'Time'
set ylabel 'Flow rate'
unset key
plot 'log' u 2:4 w l

! rm -f field-?
