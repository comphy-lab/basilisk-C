! awk '{ if ($1 == "file:") file = $2; else print $0 > file; }' < out
 
set term @PNG enhanced size 640,640

set size ratio -1
set xrange [-0.5:0.5]
set yrange [-0.5:0.5]
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

s=1.6
dry=1e-4
set multiplot layout 3,3 scale s,s
splot 'eta-0' u 1:2:($3 > dry ? $3 : 1e1000)
splot 'eta-1' u 1:2:($3 > dry ? $3 : 1e1000)
splot 'eta-2' u 1:2:($3 > dry ? $3 : 1e1000)
splot 'eta-3' u 1:2:($3 > dry ? $3 : 1e1000)
splot 'eta-4' u 1:2:($3 > dry ? $3 : 1e1000)
splot 'eta-5' u 1:2:($3 > dry ? $3 : 1e1000)
splot 'eta-6' u 1:2:($3 > dry ? $3 : 1e1000)
splot 'eta-7' u 1:2:($3 > dry ? $3 : 1e1000)
splot 'eta-8' u 1:2:($3 > dry ? $3 : 1e1000)
unset multiplot

set output 'level.png'

set multiplot layout 3,3 scale s,s
splot 'level-0'
splot 'level-1'
splot 'level-2'
splot 'level-3'
splot 'level-4'
splot 'level-5'
splot 'level-6'
splot 'level-7'
splot 'level-8'
unset multiplot

! rm -f eta-? level-?
