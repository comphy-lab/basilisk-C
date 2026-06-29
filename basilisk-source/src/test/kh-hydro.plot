unset surface
set pm3d map
set contour base
set cntrparam levels incremental -0.4,0.1,0.4
set cntrlabel onecolor
unset ytics
unset key
splot 'out' u 1:2:3 w l lc rgbcolor "black"
