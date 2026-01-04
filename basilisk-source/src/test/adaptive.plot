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
splot 'f-0'
splot 'f-1'
splot 'f-2'
splot 'f-3'
splot 'f-4'
splot 'f-5'
unset multiplot

set output 'level.png'

set multiplot layout 2,3 scale 1.6,1.6
splot 'level-0'
splot 'level-1'
splot 'level-2'
splot 'level-3'
splot 'level-4'
splot 'level-5'
unset multiplot

set term @PNG enhanced size 480,480
set colorbox
set output 'error.png'

splot 'error'

! rm -f f-? level-?
