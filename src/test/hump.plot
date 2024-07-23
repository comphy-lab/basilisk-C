! awk '{ if ($1 == "file:") file = $2; else print $0 > file; }' < out
 
set term @PNG enhanced size 800,600

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
set xrange [-0.5:1.5]
set yrange [-0.5:0.5]

s=1.7
set multiplot layout 3,2 scale s,s
splot 'eta-0'
splot 'eta-1'
splot 'eta-2'
splot 'eta-3'
splot 'eta-4'
unset multiplot

set output 'level.png'

set multiplot layout 3,2 scale s,s
splot 'level-0'
splot 'level-1'
splot 'level-2'
splot 'level-3'
splot 'level-4'
unset multiplot

! rm -f eta-? level-?
